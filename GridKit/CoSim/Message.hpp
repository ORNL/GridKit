#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace GridKit
{
  namespace CoSim
  {
    /// Disposition of the co-simulation adapter.
    enum class Disposition
    {
      /// Adapter is a server.
      SERVER,

      /// Adapter is a client.
      CLIENT,
    };

    /// Kinds of messages.
    enum class MessageType
    {
      /// Start a new simulation session.
      INIT,

      /// Advance simulation time to a target using inputs, then return outputs.
      STEP,

      /// Recalculate the interface at the current time without advancing it.
      UPDATE,

      /// Accept the interface point and write pending monitor samples.
      COMMIT,

      /// Finish the current study and release its model.
      END,

      /// Request shutting down a server after completing a batch.
      SHUTDOWN,

      /// Acknowledge receipt of the final response before closing a connection.
      CLOSE,

      /// Report an error or cancellation. Delivery is best effort.
      FAIL,
    };

    /// Signals exchanged between co-simulation participants.
    template <typename scalar_type>
    using SignalValues = std::map<std::string, scalar_type>;

    /// Metadata within a co-simulation message.
    ///
    /// This includes all information aside from signal values.
    struct MessageMetadata
    {
      /// Message kind.
      MessageType type{MessageType::INIT};

      /// Session this message belongs to.
      std::uint64_t session{0};

      /// Simulation time associated with the message.
      double time{0.0};

      /// Final simulation time, exchanged during INIT.
      double final_time{0.0};

      /// Communication interval, exchanged during INIT.
      double interval{0.0};

      /// Event times, exchanged during INIT to coordinate communication points.
      std::vector<double> event_times;

      /// Explanation of an error reported by a peer.
      std::string error;

      MessageMetadata() = default;

      MessageMetadata(
        MessageType type, std::uint64_t session, double time = 0.0
      )
        : type(type), session(session), time(time)
      {
      }
    };

    /// Full co-simulation message.
    template <typename scalar_type>
    struct MessageData : MessageMetadata
    {
      /// Signals contained within this message.
      SignalValues<scalar_type> signals;

      MessageData() = default;

      MessageData(MessageType type, std::uint64_t session, double time = 0.0)
        : MessageMetadata(type, session, time)
      {
      }
    };

    /// Callback form to apply received signals and control simulation
    /// continuity.
    template <typename scalar_type>
    using StepFunction = std::function<bool(const SignalValues<scalar_type>&)>;

    /// Callback form to extract signal data from a simulation.
    template <typename scalar_type>
    using ExtractFunction = std::function<SignalValues<scalar_type>()>;

    /// Function to send a message and wait for the next incoming message.
    ///
    /// With the client disposition, this is a request then a reply. For a
    /// server, this is a reply to the current request and the next request.
    template <typename scalar_type>
    using MessageExchangeFunction =
      std::function<MessageData<scalar_type>(const MessageData<scalar_type>&)>;

    /// Apply an incoming message and control simulation continuity.
    template <typename scalar_type>
    using MessageStepFunction = std::function<bool(
      const MessageData<scalar_type>&,
      const MessageExchangeFunction<scalar_type>&
    )>;

    /// Produce the next outgoing message, optionally exchanging intermediate
    /// messages first.
    template <typename scalar_type>
    using MessageExtractFunction = std::function<
      MessageData<scalar_type>(const MessageExchangeFunction<scalar_type>&)>;

    /// Result of a run.
    enum class RunResult
    {
      /// Run completed.
      COMPLETED,

      /// Run was stopped.
      STOPPED
    };
  } // namespace CoSim
} // namespace GridKit
