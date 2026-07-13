#pragma once

namespace Integrator
{
  /**
   * @brief Define control flow for `StepController`s to be able to control the step size of a `Rosenbrock` integrator.
   *
   */
  template <typename RealT>
  struct StepControl
  {
    /**
     * @brief Whether or not the step is accepted. A rejected step will cause the time step controller to discard
     *        the next state and re-step with the new `step_size`.
     *
     */
    bool  accept_;
    /**
     * @brief The step size the next step should take.
     *
     */
    RealT step_size_;
  };
} // namespace Integrator
