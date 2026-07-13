#include "FixedStep.hpp"

namespace Integrator
{
  /**
   * @brief Fixed step - accept every step, no matter the error, and keep the step size the same.
   *
   */
  template <typename RealT>
  StepControl<RealT> FixedStep<RealT>::nextStep([[maybe_unused]] RealT err, StepControl<RealT> prev_step, [[maybe_unused]] uint8_t method_order)
  {
    return StepControl<RealT>{
        .accept_    = true,
        .step_size_ = prev_step.step_size_,
    };
  }

  template class FixedStep<double>;
} // namespace Integrator
