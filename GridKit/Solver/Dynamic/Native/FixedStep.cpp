#include "FixedStep.hpp"

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    /**
     * @brief Fixed step - accept every step, no matter the error, and keep the step size the same.
     *
     */
    template <typename real_type>
    StepControl<real_type> FixedStep<real_type>::nextStep([[maybe_unused]] RealT err, StepControl<RealT> prev_step, [[maybe_unused]] uint8_t method_order)
    {
      return StepControl<RealT>{
          .accept_    = true,
          .step_size_ = prev_step.step_size_,
      };
    }

    template class FixedStep<double>;
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
