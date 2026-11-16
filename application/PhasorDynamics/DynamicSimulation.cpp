#include <algorithm>
#include <array>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <optional>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"
#ifdef GRIDKIT_ENABLE_ZMQ
#include <zmq.hpp>
#endif

using namespace AnalysisManager::Sundials;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int runApplication(int argc, const char* argv[])
{
  // Print summaries, such as the run time, without lowering a higher verbosity
  Log::raiseVerbosity(Log::SUMMARY);

  // Study file
  checkCommandLine(argc, "DynamicSimulation");
  auto study = parseStudyData(argv[1]);

  using EventType = SystemEvent::Type;

  for (const auto& event : study.events)
  {
    if ((event.type == EventType::FAULT_ON || event.type == EventType::FAULT_OFF)
        && event.bus == GridKit::INVALID_INDEX<index_type>)
    {
      Log::error() << "Fault event at time " << event.time
                   << " requires a bus number\n";
      return 1;
    }
  }

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);

  // Exchanges the interface bus signals with the co-simulation peers at each step
  std::optional<std::function<void(real_type)>> exchange;

#ifdef GRIDKIT_ENABLE_ZMQ
  // Signals of an interface bus and the socket to its co-simulation peer
  struct Link
  {
    bool                                                      replies{};
    zmq::socket_t                                             socket;
    std::array<scalar_type, 2>                                input{};
    std::array<index_type, 2>                                 input_index{};
    std::array<const SignalNode<scalar_type, index_type>*, 2> output{};
  };

  zmq::context_t    context;
  std::vector<Link> links(study.connections.size());

  for (std::size_t i = 0; i < links.size(); ++i)
  {
    const auto& connection = study.connections[i];
    auto&       link       = links[i];

    // The interface bus receives its input signals from the peer and sends its
    // output signals to the peer
    using BusDataT   = BusData<real_type, index_type>;
    const auto  bus  = std::ranges::find(study.model_data.bus, connection.bus, &BusDataT::bus_id);
    std::size_t port = 0;
    link.input_index.fill(GridKit::INVALID_INDEX<index_type>);
    for (const auto& signal : bus->signal_inputs)
    {
      sys.getSignalNode(signal.second)->link(&link.input[port], &link.input_index[port]);
      ++port;
    }
    port = 0;
    for (const auto& signal : bus->signal_outputs)
    {
      link.output[port] = sys.getSignalNode(signal.second);
      ++port;
    }

    // The voltage input side binds the endpoint and replies with the current
    // drawn at the voltage it receives
    const int timeout = 10000;
    link.replies      = bus->bus_type == BusDataT::BusType::SIGNAL_VOLTAGE_IN;
    if (link.replies)
    {
      link.socket = zmq::socket_t(context, zmq::socket_type::rep);
      link.socket.bind(connection.address);
    }
    else
    {
      link.socket = zmq::socket_t(context, zmq::socket_type::req);
      link.socket.connect(connection.address);
    }
    link.socket.set(zmq::sockopt::rcvtimeo, timeout);
    link.socket.set(zmq::sockopt::linger, timeout);
  }

  auto receive = [](Link& link)
  {
    if (!link.socket.recv(zmq::buffer(link.input), zmq::recv_flags::none))
    {
      throw std::runtime_error("Peer did not respond");
    }
  };

  auto send = [](Link& link)
  {
    const std::array<scalar_type, 2> values{link.output[0]->read(), link.output[1]->read()};
    link.socket.send(zmq::buffer(values), zmq::send_flags::none);
  };

  if (!links.empty())
  {
    // Send the interface values at time t and receive those of the peers. A
    // voltage input side integrates each step at the voltage its peer reached,
    // so it has nothing to receive at the final time.
    exchange = [&](real_type t)
    {
      sys.evaluateResidual();
      for (auto& link : links)
      {
        send(link);
      }
      for (auto& link : links)
      {
        if (!link.replies || t < study.tmax)
        {
          receive(link);
        }
      }
    };
  }

  // A voltage input side initializes at the voltage it receives
  for (auto& link : links)
  {
    if (link.replies)
    {
      receive(link);
    }
  }
#endif

  sys.allocate();

  if (exchange)
  {
    (*exchange)(0.0);
  }

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setMaxOrder(study.max_order);
  ida.setConsistentICType(study.consistent_ic_type);
  // Received inputs are discontinuous and the algebraic variables jump with them
  ida.setSuppressAlgebraicErrors(exchange.has_value());
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  // Initialize simulation for first run
  auto      dt_monitor = study.dt_monitor;
  real_type final_time = study.tmax;
  ida.initializeSimulation(0.0);
  for (const auto& event : study.events)
  {
    // Run to event time
    ida.runSimulation(event.time, dt_monitor, exchange);

    // Set up run for event (to start at event time)
    int ret = 0;

    switch (event.type)
    {
    case EventType::FAULT_ON:
      ret = sys.getBus(event.bus)->setFault(true, event.R, event.X);
      break;
    case EventType::FAULT_OFF:
      ret = sys.getBus(event.bus)->clearFault();
      break;
    case EventType::IN_SERVICE:
      sys.getBranch(event.device)->setInService(true);
      break;
    case EventType::OUT_OF_SERVICE:
      sys.getBranch(event.device)->setInService(false);
      break;
    }

    if (ret != 0)
    {
      Log::error() << "Fault event failed for bus " << event.bus
                   << " at time " << event.time << '\n';
      return 1;
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time);
  }

  // Run to final time
  ida.runSimulation(final_time, dt_monitor, exchange);

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  Log::summary() << "Complete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status.get();
}

int main(int argc, const char* argv[])
{
  try
  {
    return runApplication(argc, argv);
  }
  catch (const std::exception& error)
  {
    Log::error() << "DynamicSimulation failed: " << error.what() << '\n';
  }

  return 1;
}
