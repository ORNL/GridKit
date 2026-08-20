#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include "CoSim.hpp"
#include <zmq.hpp>

using namespace GridKit;
using namespace GridKit::PhasorDynamics;
using namespace AnalysisManager::Sundials;

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
  using SignalT      = typename SystemModelT::SignalT;

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
    ir_signal_->set(&ir_, &ir_idx_);
    ii_signal_->set(&ii_, &ii_idx_);
    socket_.connect("tcp://0.0.0.0:5556");
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
  }

  /**
   * @brief Send voltage to and receive current from server-side instance for a
   * single time step.
   */
  void exchange()
  {
    // 1. Send data
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(16);
    oss << CoSim::STEP << " " << vr_signal_->read() << " " << vi_signal_->read();
    // std::cout << "[CLIENT] Sending: " << oss.str() << std::endl;

    zmq::message_t s_msg{oss.str().data(), oss.str().size()};
    socket_.send(s_msg, zmq::send_flags::none);

    // 2. Receive data from DataBroker
    zmq::message_t r_msg;
    auto           recv_result = socket_.recv(r_msg, zmq::recv_flags::none);
    if (recv_result)
    {
      std::istringstream iss(r_msg.to_string());
      // std::cout << "[CLIENT] Received: " << iss.str() << std::endl;
      iss >> ir_ >> ii_;
    }
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
  // Instantiate system
  auto filepath = std::filesystem::path("ThreeBusCoSimClient.case.json");
  auto data     = parseSystemModelData(filepath);
  auto sys      = SystemModel<ScalarT, IdxT>(data);
  auto client   = CoSimClient<ScalarT, IdxT>(
      sys.getSignal(1), sys.getSignal(2), sys.getSignal(3), sys.getSignal(4));
  sys.allocate();
  client.exchange();

  // Set up simulation
  Ida<ScalarT, IdxT> ida(&sys);
  ida.setTolerance(1.0e-7, 1.0e-9);
  ida.configureSimulation();

  // TODO: Take one step at a time and exchange data between.
  //       Use step_callback for now.
  auto step_cb = [&client](auto)
  {
    client.exchange();
  };

  RealT dt = 1.0 / 4.0 / 60.0;

  // Run for 1s
  ida.initializeSimulation(0.0);
  ida.runSimulation(1.0, dt, step_cb);

  // Introduce fault and run for the next 0.1s
  sys.getBusFault(0)->setStatus(true);
  ida.initializeSimulation(1.0);
  ida.runSimulation(1.1, dt, step_cb);

  // Clear the fault and run until t = 10s.
  sys.getBusFault(0)->setStatus(false);
  ida.initializeSimulation(1.1);
  ida.runSimulation(10.0, dt, step_cb);

  return 0;
}
