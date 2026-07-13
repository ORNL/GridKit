#pragma once

#include <GridKit/Solver/Dynamic/StepController.hpp>

namespace Integrator
{
  /**
   * @brief A fixed step controller which doesn't change the step size and accepts every step.
   *        Useful if you know what time scale your simulation operates on apriori and you're
   *        using a method without an embedded error controller.
   *
   *        To set the fixed size, set the `Rosenbrock::Parameters::starting_step` parameter.
   *
   */
  template <typename RealT>
  class FixedStep : public StepController<RealT>
  {
    StepControl<RealT> nextStep(RealT err, StepControl<RealT> prev_step, uint8_t method_order) final;

    /**
     * @brief This controller does not use error estimates.
     *
     * @see `nextStep()`
     *
     */
    constexpr bool usesError() const final
    {
      return false;
    }
  };
} // namespace Integrator