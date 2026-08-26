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
  using SignalNodeT  = typename SystemModelT::SignalNodeT;
  /// Alias for Ida
  using IdaT         = Ida<ScalarT, IdxT>;

  CoSimServer() = delete;

  /**
   * @brief Construct with set of signal nodes to connect
   *
   * This also binds the tcp port to which the client is expected to connect
   *
   * @param sr_out node from which to read real component of signal to send
   * @param si_out node from which to read imaginary component of signal to send
   * @param sr_in node for communicating received real component of signal
   * @param si_in node for communicating received imaginary component of signal
   */
  CoSimServer(SystemModelT& sys, SignalNodeT* sr_out, SignalNodeT* si_out, SignalNodeT* sr_in, SignalNodeT* si_in)
    : sr_out_(sr_out),
      si_out_(si_out),
      sr_in_(sr_in),
      si_in_(si_in),
      ida_(&sys),
      ctx_{},
      socket_(ctx_, zmq::socket_type::rep)
  {
    sr_in_->link(&sr_, &sr_idx_);
    si_in_->link(&si_, &si_idx_);

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
    ida_.setTolerance(1.0e-5, 1.0e-7);
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
      if (status == CoSim::Status::FAIL)
      {
        throw std::runtime_error("Received failure message from client");
      }

      // 2. Perform action
      try
      {
        if (status == CoSim::Status::INIT)
        {
          iss >> ti_ >> tf_ >> dt_ >> nsteps_ >> sr_ >> si_;
          ida_.initializeSimulation(ti_);
        }
        else if (status == CoSim::Status::STEP)
        {
          iss >> step_ >> sr_ >> si_;
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
      auto sr = sr_out_->read();
      auto si = si_out_->read();
      oss << CoSim::Status::STEP << " " << sr << " " << si;
      Log::misc() << "SERVER: Sending \"" << oss.str() << "\"\n";

      zmq::message_t reply{oss.str().data(), oss.str().size()};
      socket_.send(reply, zmq::send_flags::none);

    } while (status != CoSim::Status::END);

    Log::summary() << "SERVER: Simulation stopped\n";
  }

private:
  /// node from which to read real component of signal to send
  SignalNodeT* sr_out_;
  /// node from which to read imaginary component of signal to send
  SignalNodeT* si_out_;
  /// node for communicating received real component of signal
  SignalNodeT* sr_in_;
  /// node for communicating received imaginary component of signal
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

  auto server = CoSimServer<ScalarT, IdxT>(sys,
                                           sys.getSignalNode(1),
                                           sys.getSignalNode(2),
                                           sys.getSignalNode(3),
                                           sys.getSignalNode(4));
  sys.allocate();

  server.start();

  return 0;
}
