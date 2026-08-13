#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

#include "CoSim.hpp"
#include <zmq.hpp>

using ScalarT = double;
using IdxT    = std::size_t;

using namespace GridKit;
using namespace GridKit::PhasorDynamics;

template <typename scalar_type, typename index_type>
class CoSimServer
{
public:
  using ScalarT      = scalar_type;
  using IdxT         = index_type;
  using SystemModelT = GridKit::PhasorDynamics::SystemModel<ScalarT, IdxT>;
  using SignalT      = typename SystemModelT::SignalT;

  CoSimServer() = delete;

  CoSimServer(SignalT* ir, SignalT* ii)
    : ir_signal_(ir),
      ii_signal_(ii),
      ctx_{},
      socket_(ctx_, zmq::socket_type::rep)
  {
    socket_.bind("tcp://0.0.0.0:5556");
  }

  void start()
  {
    CoSim::StatusT status;
    ScalarT        d1, d2;
    do
    {
      // 1. Receive data
      zmq::message_t msg;
      auto           recv_result = socket_.recv(msg, zmq::recv_flags::none);
      if (recv_result)
      {
        std::istringstream(msg.to_string()) >> status >> d1 >> d2;
        // std::cout << "[SERVER] Received: " << msg.to_string_view() << std::endl;
      }

      // 2. Respond with new data
      std::ostringstream oss;
      oss << std::scientific << std::setprecision(16);
      oss << ir_signal_->read() << " " << ii_signal_->read();
      // std::cout << "[SERVER] Sending: " << oss.str() << std::endl;

      zmq::message_t reply{oss.str().data(), oss.str().size()};
      socket_.send(reply, zmq::send_flags::none);
    } while (status == CoSim::STEP);
  }

private:
  SignalT*       ir_signal_;
  SignalT*       ii_signal_;
  zmq::context_t ctx_;
  zmq::socket_t  socket_;
};

int main()
{
  // Instantiate system
  auto data = parseSystemModelData("ThreeBusCoSimServer.case.json");
  auto sys  = SystemModel<ScalarT, IdxT>(data);
  sys.allocate();

  // Set up cosim
  CoSimServer<ScalarT, IdxT> server(sys.getSignal(1), sys.getSignal(2));
  server.start();

  return 0;
}
