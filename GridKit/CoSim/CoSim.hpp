#pragma once

#include <atomic>
#include <string>

#include <GridKit/CoSim/Message.hpp>

namespace GridKit
{
  namespace CoSim
  {
    namespace detail
    {
      /// State of a co-simulation run.
      enum class RunState : std::uint8_t
      {
        /// Simulation is not running.
        IDLE,

        /// Simulation is running.
        RUNNING,

        /// Simulation stop has been requested.
        STOP_REQUESTED
      };
    } // namespace detail

    /// Co-simulation adapter.
    template <typename scalar_type, Disposition disposition>
    class CoSim
    {
      static_assert(
        std::atomic<detail::RunState>::is_always_lock_free,
        "lock-free atomics for std::uint8_t are required for co-simulation"
      );

    public:
      CoSim(
        const std::string&           address,
        StepFunction<scalar_type>    step,
        ExtractFunction<scalar_type> extract,
        int                          timeout_ms = 30000
      );
      CoSim(
        const std::string&                  address,
        MessageStepFunction<scalar_type>    step,
        MessageExtractFunction<scalar_type> extract,
        int                                 timeout_ms = 30000
      );
      CoSim(const CoSim&)            = delete;
      CoSim& operator=(const CoSim&) = delete;

      /// Run the co-simulation.
      RunResult run();

      /// Request the co-simulation to stop.
      void stop() noexcept
      {
        auto expected = detail::RunState::RUNNING;
        state_.compare_exchange_strong(
          expected, detail::RunState::STOP_REQUESTED, std::memory_order_relaxed
        );
      }

      /// Check if the co-simulation is running.
      bool running() const noexcept
      {
        return state_.load(std::memory_order_acquire) != detail::RunState::IDLE;
      }

    private:
      /// Endpoint address of the co-simulation.
      ///
      /// If the adapter disposition is `Disposition::SERVER`, this
      /// represents the address the adapter will listen on. Otherwise,
      /// this is the address of a co-simulation server to connect to.
      std::string address_;

      /// Step function.
      MessageStepFunction<scalar_type> step_;

      /// Extraction function.
      MessageExtractFunction<scalar_type> extract_;

      /// Timeout of ZeroMQ operations in milliseconds.
      int timeout_ms_;

      /// Current state of the co-simulation adapter.
      std::atomic<detail::RunState> state_{detail::RunState::IDLE};
    };

    extern template class CoSim<double, Disposition::SERVER>;
    extern template class CoSim<double, Disposition::CLIENT>;
    extern template class CoSim<float, Disposition::SERVER>;
    extern template class CoSim<float, Disposition::CLIENT>;
  } // namespace CoSim
} // namespace GridKit
