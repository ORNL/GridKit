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
  using SignalT      = typename SystemModelT::SignalT;
  /// Alias for Ida
  using IdaT         = Ida<ScalarT, IdxT>;

  CoSimClient() = delete;

  /**
   * @brief Construct with set of signal nodes to connect
   *
   * This links the current signal nodes with variables received from server and
   * connects to the expected tcp port from which to receive.
   *
   * @param vr node from which to read real component of voltage to send
   * @param vi node from which to read imaginary component of voltage to send
   * @param ir node for communicating received real component of current
   * @param ii node for communicating received imaginary component of current
   */
  CoSimClient(SystemModelT& sys, SignalT* vr, SignalT* vi, SignalT* ir, SignalT* ii)
    : vr_signal_(vr),
      vi_signal_(vi),
      ir_signal_(ir),
      ii_signal_(ii),
      ida_(&sys),
      ctx_{},
      socket_(ctx_, zmq::socket_type::req)
  {
    ir_signal_->set(&ir_, &ir_idx_);
    ii_signal_->set(&ii_, &ii_idx_);

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
   * @brief Send voltage to and receive current from server-side instance for a
   * single time step.
   */
  void exchange(CoSim::Status stat)
  {
    // 1. Send data
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(16);
    if (stat == CoSim::Status::INIT)
    {
      oss << stat << " " << ti_ << " " << tf_ << " " << dt_ << " " << nsteps_;
    }
    else
    {
      auto vr = vr_signal_->read();
      auto vi = vi_signal_->read();
      oss << stat << " " << step_ << " " << vr << " " << vi;
    }

    Log::misc() << "CLIENT: Sending \"" << oss.str() << "\"\n";
    zmq::message_t s_msg{oss.str().data(), oss.str().size()};
    socket_.send(s_msg, zmq::send_flags::none);

    // 2. Receive data from DataBroker
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
    iss >> ir_ >> ii_;
  }

  void runSimulation(RealT ti, RealT tf, RealT dt)
  {
    ti_ = ti;
    tf_ = tf;
    dt_ = dt;

    ida_.initializeSimulation(ti_);

    nsteps_ = ida_.getStepCount(tf_, dt_);

    exchange(CoSim::Status::INIT);

    for (step_ = 1; step_ <= nsteps_; step_++)
    {
      ida_.runSimulationStep(tf_, dt_, step_, nsteps_);
      exchange(CoSim::Status::STEP);
    }
  }

private:
  /// node from which to read real component of voltage to send
  SignalT* vr_signal_;
  /// node from which to read imaginary component of voltage to send
  SignalT* vi_signal_;
  /// node for communicating received real component of current
  SignalT* ir_signal_;
  /// node for communicating received imaginary component of current
  SignalT* ii_signal_;
  /// variable for receiving current
  ScalarT  ir_{};
  /// variable for receiving current
  ScalarT  ii_{};
  /// dummy index for current signal
  IdxT     ir_idx_{GridKit::INVALID_INDEX<IdxT>};
  /// dummy index for current signal
  IdxT     ii_idx_{GridKit::INVALID_INDEX<IdxT>};

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
  auto data     = parseSystemModelData(filepath);
  auto sys      = SystemModel<ScalarT, IdxT>(data);

  auto client = CoSimClient<ScalarT, IdxT>(sys,
                                           sys.getSignal(1),
                                           sys.getSignal(2),
                                           sys.getSignal(3),
                                           sys.getSignal(4));
  sys.allocate();
  client.configure();

  RealT dt = 1.0 / 4.0 / 60.0;

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
