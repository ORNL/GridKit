#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSim.hpp"
#include <zmq.hpp>

using ScalarT = double;
using IdxT    = std::size_t;

using namespace GridKit;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Utilities;
using namespace AnalysisManager::Sundials;

using Log = GridKit::Utilities::Logger;

/**
 * @brief A simple implementation of the "server" side of a co-simulation pair
 *
 * This is a "server" in the sense that it waits for a request from the "client"
 * to initiate each step in the simulation and to trigger when to stop the
 * simulation.
 *
 * @tparam scalar_type scalar parameter type
 * @tparam index_type  integer parameter type
 */
template <typename scalar_type, typename index_type>
class CoSimServer
{
public:
  /// Type representing a scalar value
  using ScalarT      = scalar_type;
  /// Type representing an index
  using IdxT         = index_type;
  /// Type representing a real value
  using RealT        = ScalarTraits<ScalarT>::RealT;
  /// Alias for SystemModel
  using SystemModelT = SystemModel<ScalarT, IdxT>;
  /// Alias for BusBase
  using BusT         = typename SystemModelT::BusT;
  /// Alias for SignalNode
  using SignalT      = typename SystemModelT::SignalT;
  /// Alias for Ida
  using IdaT         = Ida<ScalarT, IdxT>;

  CoSimServer() = delete;

  /**
   * @brief Construct with set of signal nodes to connect
   *
   * This also binds the tcp port to which the client is expected to connect
   *
   * @param bus interfacing bus
   * @param ir node from which to read real component of current to send
   * @param ii node from which to read imaginary component of current to send
   */
  CoSimServer(SystemModelT& sys, BusT* bus)
    : bus_(bus),
      ida_(&sys),
      ctx_{},
      socket_(ctx_, zmq::socket_type::rep)
  {
    socket_.bind("tcp://0.0.0.0:5556");
  }

  /**
   * @brief Start the server
   *
   * The server will stay in this function until the end of the simulation is
   * triggered by the client.
   *
   * Each time a voltage message is received, a step is taken and the resulting
   * currents are sent as a response.
   *
   * Messages received from the client begin with a status token. Stepping will
   * continue as long as the status received is CoSim::STEP. When CoSim::END is
   * received the simulation will wrap up.
   */
  void start()
  {
    Log::summary() << "SERVER: Start simulation loop\n";
    ida_.setTolerance(1.0e-7, 1.0e-9);
    ida_.configureSimulation();

    CoSim::Status status;

    do
    {
      // 1. Receive data
      zmq::message_t r_msg;
      auto           recv_result = socket_.recv(r_msg, zmq::recv_flags::none);
      if (!recv_result)
      {
        throw std::runtime_error("Co-sim server: failed receive");
      }
      auto r_msg_str = r_msg.to_string();
      Log::misc() << "SERVER: Received \"" << r_msg_str << "\"\n";
      auto iss = std::istringstream(r_msg_str);
      iss >> status;

      // 2. Perform action
      try
      {
        if (status == CoSim::Status::INIT)
        {
          iss >> ti_ >> tf_ >> dt_ >> nsteps_;
          ida_.initializeSimulation(ti_);
        }
        else if (status == CoSim::Status::STEP)
        {
          iss >> step_ >> bus_->Vr() >> bus_->Vi();
          ida_.runSimulationStep(tf_, dt_, step_, nsteps_);
        }
      }
      catch (const std::exception& e)
      {
        // 3. Send failure
        std::ostringstream oss;
        oss << CoSim::Status::FAIL;
        zmq::message_t fail_msg{oss.str().data(), oss.str().size()};
        socket_.send(fail_msg, zmq::send_flags::none);
        throw;
      }

      // 3. Respond with new data
      std::ostringstream oss;
      oss << std::scientific << std::setprecision(16);
      oss << CoSim::Status::STEP << " " << bus_->Ir() << " " << bus_->Ii();
      Log::misc() << "SERVER: Sending \"" << oss.str() << "\"\n";

      zmq::message_t reply{oss.str().data(), oss.str().size()};
      socket_.send(reply, zmq::send_flags::none);

    } while (status != CoSim::Status::END);

    Log::summary() << "SERVER: Simulation stopped\n";
  }

private:
  /// interfacing bus
  BusT* bus_;

  /// Ida solver
  IdaT ida_;

  // Solver stepping parameters
  RealT ti_;
  RealT tf_;
  RealT dt_;
  int   step_;
  int   nsteps_;

  /// ZMQ context
  zmq::context_t ctx_;
  /// ZMQ socket
  zmq::socket_t  socket_;
};

int main(int argc, const char* argv[])
{
  Log::setVerbosity(Log::Verbosity::EVERYTHING);

  CliArgs args{{.name     = {"--case-file", "-c"},
                .required = true},

               {.name = {"--iface-bus", "-b"},
                .type = ArgType::Integer}};

  args.parseArgs(argc, argv);

  auto filepath = std::filesystem::path(args["case-file"]());
  auto busid    = args["iface-bus"].as<IdxT>();

  // Instantiate system
  auto data = parseSystemModelData(filepath);
  auto sys  = SystemModel<ScalarT, IdxT>(data);

  auto server = CoSimServer<ScalarT, IdxT>(sys, sys.getBus(busid));
  sys.allocate();

  server.start();

  return 0;
}
