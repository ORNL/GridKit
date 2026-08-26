#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSim.hpp"
#include <zmq.hpp>

using namespace GridKit;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Utilities;
using namespace AnalysisManager::Sundials;

using Log = GridKit::Utilities::Logger;

/**
 * @brief A simple implementation of the "client" side of a co-simulation pair
 *
 * This is a "client" in the sense that it initiates the simulation and triggers
 * its end.
 *
 * @tparam scalar_type scalar parameter type
 * @tparam index_type  integer parameter type
 */
template <typename scalar_type, typename index_type>
class CoSimClient
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
  /// Alias for SignalNode
  using SignalNodeT  = typename SystemModelT::SignalNodeT;
  /// Alias for Ida
  using IdaT         = Ida<ScalarT, IdxT>;

  CoSimClient() = delete;

  /**
   * @brief Construct with set of signal nodes to connect
   *
   * This links the current signal nodes with variables received from server and
   * connects to the expected tcp port from which to receive.
   *
   * @param sr_out node from which to read real component of signal to send
   * @param si_out node from which to read imaginary component of signal to send
   * @param sr_in node for communicating received real component of signal
   * @param si_in node for communicating received imaginary component of signal
   */
  CoSimClient(SystemModelT& sys, SignalNodeT* sr_out, SignalNodeT* si_out, SignalNodeT* sr_in, SignalNodeT* si_in)
    : sr_out_(sr_out),
      si_out_(si_out),
      sr_in_(sr_in),
      si_in_(si_in),
      ida_(&sys),
      ctx_{},
      socket_(ctx_, zmq::socket_type::req)
  {
    sr_in_->link(&sr_, &sr_idx_);
    si_in_->link(&si_, &si_idx_);

    socket_.connect("tcp://0.0.0.0:5556");
    Log::summary() << "CLIENT: Established connection with server\n";
  }

  /**
   * @brief Signal the end of the simulation to the server and destruct object
   */
  ~CoSimClient()
  {
    exchange(CoSim::Status::END);
    Log::summary() << "CLIENT: Ending simulation\n";
  }

  /**
   * @brief Perform initial solver setup
   */
  void configure()
  {
    ida_.setTolerance(1.0e-7, 1.0e-9);
    ida_.configureSimulation();
  }

  /**
   * @brief Send output signals to and receive input signals from server-side
   * instance for a single time step.
   */
  void exchange(CoSim::Status stat)
  {
    // 1. Send data
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(16);
    oss << stat;
    if (stat == CoSim::Status::INIT)
    {
      auto sr = sr_out_->read();
      auto si = si_out_->read();
      oss << " " << ti_ << " " << tf_ << " " << dt_ << " " << nsteps_
          << " " << sr << " " << si;
    }
    else if (stat != CoSim::Status::FAIL)
    {
      auto sr = sr_out_->read();
      auto si = si_out_->read();
      oss << " " << step_ << " " << sr << " " << si;
    }

    Log::misc() << "CLIENT: Sending \"" << oss.str() << "\"\n";
    zmq::message_t s_msg{oss.str().data(), oss.str().size()};
    socket_.send(s_msg, zmq::send_flags::none);
    if (stat == CoSim::Status::FAIL)
    {
      return;
    }

    // 2. Receive data
    zmq::message_t r_msg;
    auto           recv_result = socket_.recv(r_msg, zmq::recv_flags::none);
    if (!recv_result)
    {
      throw std::runtime_error("Co-sim client: failed receive");
    }
    auto r_msg_str = r_msg.to_string();
    Log::misc() << "CLIENT: Received \"" << r_msg_str << "\"\n";
    auto iss = std::istringstream(r_msg_str);

    CoSim::Status rstat;
    iss >> rstat;
    if (rstat == CoSim::Status::FAIL)
    {
      throw std::runtime_error("Received failure message from server");
    }
    iss >> sr_ >> si_;
  }

  /**
   * @brief Run Ida over specified time period
   * @param ti Initial time
   * @param tf Final time
   * @param dt Time step size
   */
  void runSimulation(RealT ti, RealT tf, RealT dt)
  {
    ti_ = ti;
    tf_ = tf;
    dt_ = dt;

    ida_.initializeSimulation(ti_);

    nsteps_ = ida_.getStepCount(tf_, dt_);

    exchange(CoSim::Status::INIT);

    try
    {
      for (step_ = 1; step_ <= nsteps_; step_++)
      {
        ida_.runSimulationStep(tf_, dt_, step_, nsteps_);
        exchange(CoSim::Status::STEP);
      }
    }
    catch (const std::exception& e)
    {
      exchange(CoSim::Status::FAIL);
      throw;
    }
  }

private:
  /// node from which to read real component of signal to send
  SignalNodeT* sr_out_;
  /// node from which to read imaginary component of signal to send
  SignalNodeT* si_out_;
  /// node for communicating received real component of input signal
  SignalNodeT* sr_in_;
  /// node for communicating received imaginary component of input signal
  SignalNodeT* si_in_;
  /// variable for receiving signal real component
  ScalarT      sr_{};
  /// variable for receiving signal imaginary component
  ScalarT      si_{};
  /// dummy index for signal id
  IdxT         sr_idx_{GridKit::INVALID_INDEX<IdxT>};
  /// dummy index for signal id
  IdxT         si_idx_{GridKit::INVALID_INDEX<IdxT>};

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

using ScalarT = double;
using RealT   = double;
using IdxT    = std::size_t;

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

  auto client = CoSimClient<ScalarT, IdxT>(sys,
                                           sys.getSignalNode(1),
                                           sys.getSignalNode(2),
                                           sys.getSignalNode(3),
                                           sys.getSignalNode(4));
  sys.allocate();
  client.configure();

  RealT dt = 1.0 / 16.0 / 60.0;

  client.runSimulation(0.0, 0.1, dt);

  // Run for 1s
  // client.runSimulation(0.0, 1.0, dt);

  // // Introduce fault and run for the next 0.1s
  // sys.getBusFault(0)->setStatus(true);
  // client.runSimulation(1.0, 1.1, dt);

  // // Clear the fault and run until t = 10s.
  // sys.getBusFault(0)->setStatus(false);
  // client.runSimulation(1.1, 10.0, dt);

  return 0;
}
