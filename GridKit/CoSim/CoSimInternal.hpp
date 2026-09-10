#pragma once

#include <chrono>
#include <exception>
#include <stdexcept>

#include <GridKit/CoSim/CoSim.hpp>

#include <zmq.hpp>

namespace GridKit
{
  namespace CoSim
  {
    namespace detail
    {
      /// Co-simulation was stopped.
      class Stopped : public std::exception
      {
      public:
        const char* what() const noexcept override
        {
          return "co-simulation stopped";
        }
      };

      /// Co-simulation peer failed.
      class PeerFailure : public std::runtime_error
      {
      public:
        using std::runtime_error::runtime_error;
      };

      /// Private transport implementation.
      template <typename scalar_type>
      class Channel
      {
      public:
        Channel(
          Disposition                  disposition,
          const std::string&           endpoint,
          int                          timeout_ms = 30000,
          const std::atomic<RunState>* state      = nullptr
        );
        Channel(const Channel&)            = delete;
        Channel& operator=(const Channel&) = delete;

        /// Send a co-simulation message.
        void send(const MessageData<scalar_type>& message);

        /// Receive a co-simulation message.
        MessageData<scalar_type> receive();

        /// Exchange co-simulation messages with the peer.
        MessageData<scalar_type> exchange(
          const MessageData<scalar_type>& request
        );

        /// Notify the peer of failure.
        void fail(const MessageData<scalar_type>& message) noexcept;

        /// Acknowledge the last peer response prior to closing the connection.
        void closePeer(const MessageMetadata& last) noexcept;

        /// Send the final response and briefly wait for the peer to consume it.
        void finish(const MessageData<scalar_type>& reply);

        /// Bound cleanup without discarding a queued failure reply.
        void finishCancellation() noexcept;

      private:
        /// Scope guard to update the timeout to a shorter teardown timeout and
        /// restore it when the scope is exited.
        class Teardown;

        /// Wait for a ZeroMQ event.
        void wait(
          short                                 events,
          const char*                           operation,
          std::chrono::steady_clock::time_point deadline
        );

        /// Check if the simulation has been stopped.
        void checkStop() const;

        /// Wait for a close reply.
        void awaitClose();

        /// ZeroMQ context.
        zmq::context_t context_;

        /// ZeroMQ socket.
        zmq::socket_t socket_;

        /// Timeout of ZeroMQ operations in milliseconds.
        int timeout_ms_;

        /// Reference to the co-simulation state.
        const std::atomic<RunState>* state_;

        /// Whether or not we are currently able to send a message.
        bool can_send_;
      };

      extern template class Channel<double>;
      extern template class Channel<float>;
    } // namespace detail
  } // namespace CoSim
} // namespace GridKit
