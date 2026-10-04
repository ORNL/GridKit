#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSim.hpp"
#include <zmq.hpp>

using namespace GridKit;
using namespace GridKit::PhasorDynamics;
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
  /// Alias for SystemModel
  using SystemModelT = SystemModel<ScalarT, IdxT>;
  /// Alias for SignalNode
  using SignalT      = typename SystemModelT::SignalNodeT;

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
  CoSimClient(SignalT* vr, SignalT* vi, SignalT* ir, SignalT* ii)
    : vr_signal_(vr),
      vi_signal_(vi),
      ir_signal_(ir),
      ii_signal_(ii),
      ctx_{},
      socket_(ctx_, zmq::socket_type::req)
  {
    ir_signal_->link(&ir_, &ir_idx_);
    ii_signal_->link(&ii_, &ii_idx_);
    socket_.connect("tcp://0.0.0.0:5556");
    Log::summary() << "CLIENT: Established connection with server\n";
  }

  /**
   * @brief Signal the end of the simulation to the server and destruct object
   */
  ~CoSimClient()
  {
    std::ostringstream oss;
    oss << CoSim::END << " " << vr_signal_->read() << " " << vi_signal_->read();
    zmq::message_t s_msg{oss.str().data(), oss.str().size()};
    socket_.send(s_msg, zmq::send_flags::none);

    zmq::message_t r_msg;
    auto           recv_result = socket_.recv(r_msg, zmq::recv_flags::none);
    if (recv_result)
    {
    }
    Log::summary() << "CLIENT: Ending simulation\n";
  }

  /**
   * @brief Send voltage to and receive current from server-side instance for a
   * single time step.
   *
   * @return true if either received current changed.
   */
  bool exchange()
  {
    // 1. Send data
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(16);
    oss << CoSim::STEP << " " << vr_signal_->read() << " " << vi_signal_->read();
    Log::misc() << "CLIENT: Sending \"" << oss.str() << "\"\n";

    zmq::message_t s_msg{oss.str().data(), oss.str().size()};
    socket_.send(s_msg, zmq::send_flags::none);

    // 2. Receive data from DataBroker
    zmq::message_t r_msg;
    auto           recv_result = socket_.recv(r_msg, zmq::recv_flags::none);
    if (recv_result)
    {
      std::istringstream iss(r_msg.to_string());
      Log::misc() << "CLIENT: Received \"" << iss.str() << "\"\n";
      const ScalarT ir_old = ir_;
      const ScalarT ii_old = ii_;
      iss >> ir_ >> ii_;
      return ir_ != ir_old || ii_ != ii_old;
    }
    return false;
  }

private:
  /// node from which to read real component of voltage to send
  SignalT*       vr_signal_;
  /// node from which to read imaginary component of voltage to send
  SignalT*       vi_signal_;
  /// node for communicating received real component of current
  SignalT*       ir_signal_;
  /// node for communicating received imaginary component of current
  SignalT*       ii_signal_;
  /// variable for receiving current
  ScalarT        ir_{};
  /// variable for receiving current
  ScalarT        ii_{};
  /// dummy index for current signal
  IdxT           ir_idx_{GridKit::INVALID_INDEX<IdxT>};
  /// dummy index for current signal
  IdxT           ii_idx_{GridKit::INVALID_INDEX<IdxT>};
  /// ZMQ context
  zmq::context_t ctx_;
  /// ZMQ socket
  zmq::socket_t  socket_;
};

using ScalarT = double;
using RealT   = double;
using IdxT    = std::size_t;

int main()
{
  Log::raiseVerbosity(Log::Verbosity::SUMMARY);

  // Instantiate system
  auto filepath = std::filesystem::path("ThreeBusCoSimClient.case.json");
  auto data     = parseSystemModelData(filepath);
  auto sys      = SystemModel<ScalarT, IdxT>(data);
  auto client   = CoSimClient<ScalarT, IdxT>(sys.getSignalNode(1),
                                           sys.getSignalNode(2),
                                           sys.getSignalNode(3),
                                           sys.getSignalNode(4));
  sys.allocate();

  // Set up simulation
  Ida<ScalarT, IdxT> ida(&sys);
  ida.setTolerance(1.0e-7, 1.0e-9);
  ida.configureSimulation();

  client.exchange();

  // Hold received inputs fixed between communication times.
  auto run_interval = [&ida, &client](RealT tf, IdxT nsteps)
  {
    const RealT t0 = ida.getInitialTime();
    const RealT dt = (tf - t0) / static_cast<RealT>(nsteps);

    for (IdxT step = 1; step <= nsteps; ++step)
    {
      const RealT t = step == nsteps ? tf : std::fma(static_cast<RealT>(step), dt, t0);
      ida.runSimulation(t);

      const bool changed = client.exchange();
      if (changed && step < nsteps)
      {
        ida.initializeSimulation(t);
      }
    }
    // The caller handles fault changes and restarts at interval endpoints.
  };

  // Communicate at 240 Hz until the first fault.
  ida.initializeSimulation(0.0);
  run_interval(1.0, 240);

  // Introduce fault and run for the next 0.1s
  sys.getBusFault(0)->setStatus(true);
  ida.initializeSimulation(1.0);
  run_interval(1.1, 24);

  // Clear the fault and run until t = 10s.
  sys.getBusFault(0)->setStatus(false);
  ida.initializeSimulation(1.1);
  run_interval(10.0, 2136);

  return 0;
}
