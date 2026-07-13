#pragma once

#include <GridKit/Solver/Dynamic/StepController.hpp>

namespace Integrator
{

  /**
   * @brief A simple textbook adaptive `StepController` which seeks to meet a relative and absolute tolerance
   *        based on an error estimate.
   *
   */
  template <typename RealT>
  class AdaptiveStep : public StepController<RealT>
  {
    /**
     * @brief Parameters for the step controller.
     *
     */
    struct Parameters
    {
      /**
       * @brief The minimum multiple by which the step size can be multiplied to obtain the new step size.
       *        Increasing this can allow the integrator to be slightly more conservative in selecting the step size
       *        - decreasing the number of steps taken but increasing the risk of failing the next step.
       *
       * @note Should be between 0 and 1.
       *
       */
      RealT fac_min_   = 0.2;
      /**
       * @brief The maximum multiple by which the step size can be multiplied to obtain the new step size.
       *        Decreasing this will make the integrator more conservative in selecting the step size -
       *        increasing the number of steps taken but decreasing the risk of failing the next step.
       *
       * @note Should be greater than 1.
       *
       */
      RealT fac_max_   = 5.0;
      /**
       * @brief A "fudge factor" introduced to decrease risk of failing a step. The larger the fudge factor,
       *        the more likely steps will fail, but fewer steps will be taken.
       *
       * @note Should be between 0 and 1.
       *
       */
      RealT fac_scale_ = 0.9;
    } params_;

  public:
    AdaptiveStep(const Parameters& params)
      : params_(params)
    {
    }

    StepControl<RealT> nextStep(RealT err, StepControl<RealT> prev_step, uint8_t method_order) final;

    /**
     * @brief This controller uses error estimates.
     *
     * @see `nextStep()`
     *
     */
    constexpr bool usesError() const final
    {
      return true;
    }
  };
} // namespace Integrator
