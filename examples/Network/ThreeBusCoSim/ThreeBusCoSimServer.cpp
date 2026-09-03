#include <iostream>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CoSim.hpp"
#include <zmq.hpp>

using ScalarT = double;
using IdxT    = std::size_t;

using namespace GridKit;
using namespace GridKit::PhasorDynamics;

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
  /// Alias for SystemModel
  using SystemModelT = SystemModel<ScalarT, IdxT>;
  /// Alias for SignalNode
  using SignalT      = typename SystemModelT::SignalT;

  CoSimServer() = delete;

  /**
   * @brief Construct with set of signal nodes to connect
   *
   * This also binds the tcp port to which the client is expected to connect
   *
   * @param ir node from which to read real component of current to send
   * @param ii node from which to read imaginary component of current to send
   */
  CoSimServer(SignalT* ir, SignalT* ii)
    : ir_signal_(ir),
      ii_signal_(ii),
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
   *
   * @note For this initial implementation, no solver step is taken. Constant
   * current values are sent in response.
   */
  void start()
  {
    Log::summary() << "SERVER: Start simulation loop\n";
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
        Log::misc() << "SERVER: Received \"" << msg.to_string_view() << "\"\n";
      }

      // 2. Respond with new data
      std::ostringstream oss;
      oss << std::scientific << std::setprecision(16);
      oss << ir_signal_->read() << " " << ii_signal_->read();
      Log::misc() << "SERVER: Sending \"" << oss.str() << "\"\n";

      zmq::message_t reply{oss.str().data(), oss.str().size()};
      socket_.send(reply, zmq::send_flags::none);
    } while (status == CoSim::STEP);

    Log::summary() << "SERVER: Simulation stopped\n";
  }

private:
  /// node from which to read real component of current to send
  SignalT*       ir_signal_;
  /// node from which to read imaginary component of current to send
  SignalT*       ii_signal_;
  /// ZMQ context
  zmq::context_t ctx_;
  /// ZMQ socket
  zmq::socket_t  socket_;
};

int main()
{
  Log::setVerbosity(Log::Verbosity::SUMMARY);

  // Instantiate system
  auto filepath = std::filesystem::path("ThreeBusCoSimServer.case.json");
  auto data     = parseSystemModelData(filepath);
  auto sys      = SystemModel<ScalarT, IdxT>(data);
  sys.allocate();

  // Set up cosim
  CoSimServer<ScalarT, IdxT> server(sys.getSignal(1), sys.getSignal(2));
  server.start();

  return 0;
}
