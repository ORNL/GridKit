#include <chrono>
#include <cstdint>
#include <filesystem>
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/CoSim/CoSimInternal.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    using namespace GridKit::CoSim;
    using nlohmann::json;

    template <typename scalar_type, typename index_type>
    class CoSimTests
    {
      using Channel     = detail::Channel<scalar_type>;
      using Message     = MessageData<scalar_type>;
      using Signals     = SignalValues<scalar_type>;
      using PeerFailure = detail::PeerFailure;

      /// Retrieve an endpoint to use for the co-simulation.
      static std::string endpoint()
      {
        static unsigned serial = 0;

        const auto timestamp =
          std::chrono::steady_clock::now().time_since_epoch().count();
        const auto name =
          std::format("gridkit-unit-{}-{}", timestamp, ++serial);

        return std::format(
          "ipc://{}", (std::filesystem::temp_directory_path() / name).string()
        );
      }

    public:
      TestOutcome sessions()
      {
        TestStatus outcome;

        const auto address = endpoint();
        auto       server  = std::async(
          std::launch::async,
          [&]() -> auto
          {
            bool    outcome = true;
            Channel channel(Disposition::SERVER, address, 5000);

            for (std::uint64_t session = 1; session <= 2; session++)
            {
              auto init = channel.receive();

              outcome &= init.type == MessageType::INIT;
              outcome &= init.session == session;
              outcome &= init.signals.at("voltage real") == 1.03;

              init.signals     = {{"current from generator", 2.5}};
              init.event_times = {0.025, 0.035};
              channel.send(init);

              auto step = channel.receive();

              outcome &= step.time == 0.05;
              outcome &= step.type == MessageType::STEP;

              channel.send(step);

              auto end = channel.receive();

              outcome &= end.type == MessageType::END;

              channel.send(end);
            }

            const auto stop = channel.receive();

            outcome &= stop.type == MessageType::SHUTDOWN;

            channel.finish(stop);

            return outcome;
          }
        );

        Channel client(Disposition::CLIENT, address, 5000);

        for (std::uint64_t session = 1; session <= 2; session++)
        {
          Message init(MessageType::INIT, session);
          init.final_time = 0.05;
          init.interval   = 0.05;
          init.signals    = {{"voltage real", 1.03}};

          auto reply = client.exchange(init);

          outcome *= reply.signals.at("current from generator") == 2.5;
          outcome *= reply.event_times.size() == 2;

          Message step(MessageType::STEP, session);
          step.time = 0.05;

          client.exchange(step);

          Message end(MessageType::END, session);

          client.exchange(end);
        }

        Message stop(MessageType::SHUTDOWN, 2);

        client.exchange(stop);
        client.closePeer(stop);

        outcome *= server.get();

        return outcome.report(__func__);
      }

      TestOutcome rejectedReplies()
      {
        TestStatus outcome;

        for (int kind = 0; kind < 4; kind++)
        {
          const auto address = endpoint();
          auto       server  = std::async(
            std::launch::async,
            [&]
            {
              Channel channel(Disposition::SERVER, address, 5000);
              auto    message = channel.receive();

              switch (kind)
              {
              case 0:
                message.type  = MessageType::FAIL;
                message.error = "solver failure";
                break;
              case 1:
                message.session++;
                break;
              case 2:
                message.time = 0.1;
                break;
              case 3:
                message.type = MessageType::END;
                break;
              }

              channel.finish(message);
            }
          );

          Channel client(Disposition::CLIENT, address, 5000);

          if (kind == 0)
          {
            outcome *= throws<PeerFailure>(
              [&]
              { client.exchange(Message()); }
            );
          }
          else
          {
            outcome *=
              throws<std::runtime_error>(
                [&]
                { client.exchange(Message()); }
              );
            client.closePeer(Message());
          }

          server.get();
        }

        return outcome.report(__func__);
      }

      TestOutcome malformedMessages()
      {
        TestStatus outcome;

        for (const std::string wire :
             {"{",
              "{}",
              R"({"type":"unknown","session":1,"time":0,"signals":{}})",
              R"({"type":"step","session":-1,"time":0,"signals":{}})",
              R"({"type":"step","session":1,"time":0,"signals":{"v":null}})"})
        {
          const auto address = endpoint();
          auto       server  = std::async(
            std::launch::async,
            [&]() -> auto
            {
              bool    outcome = true;
              Channel channel(Disposition::SERVER, address, 5000);

              try
              {
                channel.receive();
              }
              catch (const json::exception& _)
              {
              }
              catch (const std::runtime_error& _)
              {
              }
              catch (...)
              {
                throw;
              }

              Message failure(MessageType::FAIL, 0);
              failure.error = "malformed request";

              channel.finish(failure);

              return outcome;
            }
          );

          zmq::context_t context;
          zmq::socket_t  socket(context, zmq::socket_type::req);
          socket.set(zmq::sockopt::linger, 0);
          socket.set(zmq::sockopt::rcvtimeo, 5000);
          socket.connect(address);

          outcome *=
            socket.send(zmq::buffer(wire), zmq::send_flags::none).has_value();

          zmq::message_t reply;

          auto ok = socket.recv(reply, zmq::recv_flags::none).has_value();
          if (!ok)
          {
            outcome *= false;
          }
          else
          {
            outcome *=
              reply.to_string().find("malformed request") != std::string::npos;

            const std::string close =
              R"({"type":"close","session":0,"time":0,"signals":{}})";

            outcome *= socket.send(zmq::buffer(close), zmq::send_flags::none)
                         .has_value();

            socket.set(zmq::sockopt::rcvtimeo, 100);

            static_cast<void>(socket.recv(reply, zmq::recv_flags::none));
          }

          outcome *= server.get();
        }

        return outcome.report(__func__);
      }

      TestOutcome boundedWaits()
      {
        TestStatus outcome;

        Channel client(Disposition::CLIENT, endpoint(), 100);

        outcome *= throws<std::runtime_error>(
          [&]
          { client.exchange(Message()); }
        );

        Message invalid;
        invalid.signals = {{"bad", std::numeric_limits<double>::infinity()}};

        Channel other(Disposition::CLIENT, endpoint(), 100);

        outcome *= throws<std::runtime_error>(
          [&]
          { other.send(invalid); }
        );

        outcome *= throws<std::invalid_argument>(
          [&]
          { Channel invalid_channel(Disposition::CLIENT, endpoint(), 0); }
        );

        return outcome.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
