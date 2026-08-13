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

template <typename scalar_type, typename index_type>
class CoSimClient
{
public:
  using ScalarT      = scalar_type;
  using IdxT         = index_type;
  using SystemModelT = SystemModel<ScalarT, IdxT>;
  using SignalT      = typename SystemModelT::SignalT;

  CoSimClient() = delete;

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
  SignalT*       vr_signal_;
  SignalT*       vi_signal_;
  SignalT*       ir_signal_;
  SignalT*       ii_signal_;
  ScalarT        ir_{};
  ScalarT        ii_{};
  IdxT           ir_idx_{GridKit::INVALID_INDEX<IdxT>};
  IdxT           ii_idx_{GridKit::INVALID_INDEX<IdxT>};
  zmq::context_t ctx_;
  zmq::socket_t  socket_;
};

using ScalarT = double;
using RealT   = double;
using IdxT    = std::size_t;

int main()
{
  // Instantiate system
  auto                       data = parseSystemModelData("ThreeBusCoSimClient.case.json");
  auto                       sys  = SystemModel<ScalarT, IdxT>(data);
  CoSimClient<ScalarT, IdxT> client(
      sys.getSignal(1), sys.getSignal(2), sys.getSignal(3), sys.getSignal(4));
  sys.allocate();
  client.exchange();

  // Set up simulation
  Ida<ScalarT, IdxT> ida(&sys);
  ida.setTolerance(1.0e-7, 1.0e-9);
  ida.configureSimulation();

  // TODO: Take one step at a time and exchange data between.
  //       Just use step_callback
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
