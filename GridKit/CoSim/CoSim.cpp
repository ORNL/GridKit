#include <algorithm>
#include <chrono>
#include <cmath>
#include <format>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/CoSim/CoSimInternal.hpp>

namespace GridKit
{
  namespace CoSim
  {
    namespace
    {
      using json = nlohmann::json;

      /// Convert a message type to a string.
      ///
      /// @throw runtime_error if the message type was unknown
      std::string typeName(MessageType type)
      {
        switch (type)
        {
        case MessageType::INIT:
          return "init";
        case MessageType::STEP:
          return "step";
        case MessageType::UPDATE:
          return "update";
        case MessageType::COMMIT:
          return "commit";
        case MessageType::END:
          return "end";
        case MessageType::SHUTDOWN:
          return "shutdown";
        case MessageType::CLOSE:
          return "close";
        case MessageType::FAIL:
          return "fail";
        }

        throw std::runtime_error("invalid co-simulation message type");
      }

      /// Convert a string to a message type.
      ///
      /// @throw runtime_error if the message type was unknown
      MessageType parseType(const std::string& name)
      {
        for (auto type :
             {MessageType::INIT,
              MessageType::STEP,
              MessageType::UPDATE,
              MessageType::COMMIT,
              MessageType::END,
              MessageType::SHUTDOWN,
              MessageType::CLOSE,
              MessageType::FAIL})
        {
          if (typeName(type) == name)
          {
            return type;
          }
        }

        throw std::runtime_error("unknown co-simulation message type: " + name);
      }

      /// Validate a co-simulation message.
      ///
      /// @throws runtime_error if the co simulation signals are invalid or if a
      ///                       time value was not finite
      template <typename scalar_type>
      void validate(const MessageData<scalar_type>& message)
      {
        if (!std::isfinite(message.time) || !std::isfinite(message.final_time)
            || !std::isfinite(message.interval))
        {
          throw std::runtime_error("non-finite co-simulation time");
        }

        for (const auto& [name, value] : message.signals)
        {
          if (name.empty() || !std::isfinite(value))
          {
            throw std::runtime_error("invalid co-simulation signal: " + name);
          }
        }

        for (auto time : message.event_times)
        {
          if (!std::isfinite(time))
          {
            throw std::runtime_error("non-finite co-simulation event time");
          }
        }
      }
    } // namespace

    template <typename scalar_type>
    detail::Channel<scalar_type>::Channel(
      Disposition                  disposition,
      const std::string&           endpoint,
      int                          timeout_ms,
      const std::atomic<RunState>* state
    )
      : context_(),
        socket_(
          context_,
          disposition == Disposition::SERVER ? zmq::socket_type::rep
                                             : zmq::socket_type::req
        ),
        timeout_ms_(timeout_ms),
        state_(state),
        can_send_(disposition == Disposition::CLIENT)
    {
      if (timeout_ms <= 0 || endpoint.empty())
      {
        throw std::invalid_argument(
          "co-simulation requires an endpoint and a positive timeout"
        );
      }

      // bound socket cleanup
      socket_.set(zmq::sockopt::linger, std::min(timeout_ms, 1000));
      socket_.set(zmq::sockopt::immediate, 1);

      if (disposition == Disposition::SERVER)
      {
        socket_.bind(endpoint);
      }
      else
      {
        socket_.connect(endpoint);
      }
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::checkStop() const
    {
      if (state_
          && state_->load(std::memory_order_relaxed)
               == RunState::STOP_REQUESTED)
      {
        throw Stopped{};
      }
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::wait(
      short                                 events,
      const char*                           operation,
      std::chrono::steady_clock::time_point deadline
    )
    {
      using Clock = std::chrono::steady_clock;

      while (true)
      {
        checkStop();
        const auto remaining =
          std::chrono::duration_cast<std::chrono::milliseconds>(
            deadline - Clock::now()
          );
        if (remaining.count() <= 0)
        {
          throw std::runtime_error(
            std::format("Co-simulation {} timed out", operation)
          );
        }

        zmq::pollitem_t item{socket_.handle(), 0, events, 0};
        try
        {
          zmq::poll(
            &item, 1, std::min(remaining, std::chrono::milliseconds(25))
          );
        }
        catch (const zmq::error_t& error)
        {
          if (error.num() == EINTR)
          {
            continue;
          }
          throw;
        }
        checkStop();
        if (item.revents & events)
        {
          return;
        }
      }
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::send(
      const MessageData<scalar_type>& message
    )
    {
      checkStop();
      if (!can_send_)
      {
        throw std::logic_error("Co-simulation socket is waiting for a reply");
      }

      validate(message);
      json wire = {
        {"type", typeName(message.type)},
        {"session", message.session},
        {"time", message.time},
        {"signals", message.signals}
      };

      if (message.type == MessageType::INIT)
      {
        wire["final_time"]  = message.final_time;
        wire["interval"]    = message.interval;
        wire["event_times"] = message.event_times;
      }

      if (message.type == MessageType::FAIL)
      {
        wire["error"] = message.error;
      }

      const auto data     = wire.dump();
      const auto deadline = std::chrono::steady_clock::now()
                            + std::chrono::milliseconds(timeout_ms_);

      do
      {
        wait(ZMQ_POLLOUT, "send", deadline);
      } while (!socket_.send(zmq::buffer(data), zmq::send_flags::dontwait));

      can_send_ = false;
    }

    template <typename scalar_type>
    MessageData<scalar_type> detail::Channel<scalar_type>::receive()
    {
      checkStop();

      if (can_send_)
      {
        throw std::logic_error(
          "Co-simulation socket must send before receiving"
        );
      }

      zmq::message_t raw;
      const auto     deadline = std::chrono::steady_clock::now()
                            + std::chrono::milliseconds(timeout_ms_);
      do
      {
        wait(ZMQ_POLLIN, "receive", deadline);
      } while (!socket_.recv(raw, zmq::recv_flags::dontwait));

      can_send_       = true;
      const auto wire = json::parse(raw.to_string_view());

      MessageData<scalar_type> message;
      message.type = parseType(wire.at("type").get<std::string>());
      if (!wire.at("session").is_number_unsigned())
      {
        throw std::runtime_error("Invalid co-simulation session ID");
      }

      message.session = wire.at("session").get<std::uint64_t>();
      message.time    = wire.at("time").get<double>();
      message.signals = wire.at("signals").get<SignalValues<scalar_type>>();

      if (message.type == MessageType::INIT)
      {
        message.final_time  = wire.at("final_time").get<double>();
        message.interval    = wire.at("interval").get<double>();
        message.event_times = wire.at("event_times").get<std::vector<double>>();
      }

      if (message.type == MessageType::FAIL)
      {
        message.error = wire.at("error").get<std::string>();
      }

      validate(message);
      return message;
    }

    template <typename scalar_type>
    MessageData<scalar_type> detail::Channel<scalar_type>::exchange(
      const MessageData<scalar_type>& request
    )
    {
      send(request);
      auto reply = receive();
      if (reply.type == MessageType::FAIL)
      {
        closePeer(reply);
        throw PeerFailure("co-simulation peer failed: " + reply.error);
      }

      const bool final_signals =
        request.type == MessageType::STEP && reply.type == MessageType::END;

      if ((!final_signals && reply.type != request.type)
          || reply.session != request.session || reply.time != request.time)
      {
        throw std::runtime_error(
          "co-simulation acknowledgement has the wrong type, session or time"
        );
      }

      return reply;
    }

    template <typename scalar_type>
    class detail::Channel<scalar_type>::Teardown
    {
      Channel&                     channel_;
      const std::atomic<RunState>* state_;
      int                          timeout_ms_;

    public:
      explicit Teardown(Channel& channel) noexcept
        : channel_(channel),
          state_(channel.state_),
          timeout_ms_(channel.timeout_ms_)
      {
        channel_.state_      = nullptr;
        channel_.timeout_ms_ = 100;
      }

      ~Teardown()
      {
        channel_.state_      = state_;
        channel_.timeout_ms_ = timeout_ms_;
      }
    };

    template <typename scalar_type>
    void detail::Channel<scalar_type>::finishCancellation() noexcept
    {
      try
      {
        socket_.set(zmq::sockopt::linger, 100);
      }
      catch (...)
      {
      }
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::fail(
      const MessageData<scalar_type>& message
    ) noexcept
    {
      Teardown teardown(*this);

      try
      {
        if (can_send_)
        {
          send(message);
          if (socket_.get(zmq::sockopt::type) == ZMQ_REQ)
          {
            const auto reply = receive();
            closePeer(reply);
          }
          else
          {
            awaitClose();
          }
        }
      }
      catch (...)
      {
      }
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::closePeer(
      const MessageMetadata& last
    ) noexcept
    {
      Teardown teardown(*this);

      try
      {
        MessageData<scalar_type> close;
        close.type    = MessageType::CLOSE;
        close.session = last.session;
        close.time    = last.time;
        send(close);
        receive();
      }
      catch (...)
      {
      }
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::awaitClose()
    {
      auto close = receive();
      if (close.type != MessageType::CLOSE)
      {
        close.type  = MessageType::FAIL;
        close.error = "Co-simulation already completed";
        close.signals.clear();
      }
      send(close);
    }

    template <typename scalar_type>
    void detail::Channel<scalar_type>::finish(
      const MessageData<scalar_type>& reply
    )
    {
      send(reply);
      Teardown teardown(*this);

      try
      {
        awaitClose();
      }
      catch (...)
      {
      }
    }

    template <typename scalar_type, Disposition disposition>
    CoSim<scalar_type, disposition>::CoSim(
      const std::string&                  address,
      MessageStepFunction<scalar_type>    step,
      MessageExtractFunction<scalar_type> extract,
      int                                 timeout_ms
    )
      : address_(address),
        step_(std::move(step)),
        extract_(std::move(extract)),
        timeout_ms_(timeout_ms)
    {
      if (address.empty() || timeout_ms <= 0 || !step_ || !extract_)
      {
        throw std::invalid_argument(
          "co-simulation requires an endpoint, callbacks and a positive timeout"
        );
      }
    }

    template <typename scalar_type, Disposition disposition>
    CoSim<scalar_type, disposition>::CoSim(
      const std::string&           address,
      StepFunction<scalar_type>    step,
      ExtractFunction<scalar_type> extract,
      int                          timeout_ms
    )
      : address_(address), timeout_ms_(timeout_ms)
    {
      if (address.empty() || timeout_ms <= 0 || !step || !extract)
      {
        throw std::invalid_argument(
          "co-simulation requires an endpoint, callbacks and a positive timeout"
        );
      }

      step_ = [step = std::move(step)](const auto& message, const auto&)
      { return step(message.signals); };
      extract_ = [extract = std::move(extract)](const auto&)
      {
        MessageData<scalar_type> message;
        message.type    = MessageType::STEP;
        message.signals = extract();
        return message;
      };
    }

    template <typename scalar_type, Disposition disposition>
    RunResult CoSim<scalar_type, disposition>::run()
    {
      auto expected = detail::RunState::IDLE;
      if (!state_.compare_exchange_strong(
            expected,
            detail::RunState::RUNNING,
            std::memory_order_acquire,
            std::memory_order_relaxed
          ))
      {
        throw std::logic_error(
          "attempted to start a co-simulation while it was already running"
        );
      }

      // ensure the state is reset on scope exit
      struct Guard
      {
        std::atomic<detail::RunState>& state;

        ~Guard()
        {
          state.store(detail::RunState::IDLE, std::memory_order_release);
        }
      } guard{state_};

      detail::Channel<scalar_type> channel(
        disposition, address_, timeout_ms_, &state_
      );
      MessageMetadata current;

      struct Completed
      {
      };

      try
      {
        MessageMetadata last_reply;
        bool            peer_finished = false;
        const auto      receive       = [&]
        {
          auto message = channel.receive();
          current      = message;

          if (message.type == MessageType::FAIL)
          {
            throw std::runtime_error(
              "co-simulation client failed: " + message.error
            );
          }

          if (message.type == MessageType::CLOSE)
          {
            channel.send(message);
            throw Completed{};
          }
          return message;
        };
        const MessageExchangeFunction<scalar_type> exchange =
          [&](const auto& message)
        {
          if constexpr (disposition == Disposition::SERVER)
          {
            if (current.type == MessageType::STEP
                && message.type == MessageType::END)
            {
              channel.finish(message);
              throw Completed{};
            }
            channel.send(message);
            return receive();
          }
          else
          {
            if (peer_finished)
            {
              throw std::runtime_error(
                "co-simulation peer ended before the callback completed"
              );
            }
            current       = message;
            auto reply    = channel.exchange(message);
            last_reply    = reply;
            peer_finished = message.type == MessageType::STEP
                            && reply.type == MessageType::END;
            return reply;
          }
        };
        while (true)
        {
          if (state_.load(std::memory_order_relaxed)
              == detail::RunState::STOP_REQUESTED)
          {
            throw detail::Stopped{};
          }

          if constexpr (disposition == Disposition::SERVER)
          {
            const auto incoming = receive();
            const bool proceed  = step_(incoming, exchange);

            if (state_.load(std::memory_order_relaxed)
                == detail::RunState::STOP_REQUESTED)
            {
              throw detail::Stopped{};
            }

            auto reply = extract_(exchange);

            // terminate the step exchange
            if (!proceed && current.type == MessageType::STEP
                && reply.type == MessageType::STEP)
            {
              reply.type = MessageType::END;
            }

            if (!proceed)
            {
              channel.finish(reply);
              return RunResult::COMPLETED;
            }

            channel.send(reply);
          }
          else
          {
            auto outgoing = extract_(exchange);
            bool proceed  = false;
            if (!peer_finished)
            {
              const auto incoming = exchange(outgoing);
              proceed             = step_(incoming, exchange);
            }
            if (!proceed || peer_finished)
            {
              channel.closePeer(last_reply);
              return state_.load(std::memory_order_relaxed)
                         == detail::RunState::STOP_REQUESTED
                       ? RunResult::STOPPED
                       : RunResult::COMPLETED;
            }
          }
        }
      }
      catch (const Completed&)
      {
        return RunResult::COMPLETED;
      }
      catch (const detail::Stopped&)
      {
        MessageData<scalar_type> failure;
        failure.type    = MessageType::FAIL;
        failure.session = current.session;
        failure.time    = std::isfinite(current.time) ? current.time : 0.0;
        failure.error   = "co-simulation stopped";
        channel.fail(failure);
        channel.finishCancellation();
        return RunResult::STOPPED;
      }
      catch (const detail::PeerFailure&)
      {
        throw;
      }
      catch (const std::exception& error)
      {
        MessageData<scalar_type> failure;
        failure.type    = MessageType::FAIL;
        failure.session = current.session;
        failure.time    = std::isfinite(current.time) ? current.time : 0.0;
        failure.error   = error.what();
        channel.fail(failure);
        throw;
      }
      catch (...)
      {
        MessageData<scalar_type> failure;
        failure.type    = MessageType::FAIL;
        failure.session = current.session;
        failure.time    = std::isfinite(current.time) ? current.time : 0.0;
        failure.error   = "non-standard exception in co-simulation callback";
        channel.fail(failure);
        throw;
      }
    }

    template class CoSim<double, Disposition::SERVER>;
    template class CoSim<double, Disposition::CLIENT>;
    template class CoSim<float, Disposition::SERVER>;
    template class CoSim<float, Disposition::CLIENT>;
    template class detail::Channel<double>;
    template class detail::Channel<float>;
  } // namespace CoSim
} // namespace GridKit
