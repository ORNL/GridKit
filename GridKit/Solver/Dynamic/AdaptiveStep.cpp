#include "AdaptiveStep.hpp"

#include <cmath>

#include <GridKit/Constants.hpp>

namespace Integrator
{
  /**
   * @brief Standard textbook adaptive controller. Accept if `err <= 1` and use
   *
   * \f[h_{new} = h * \min \left\{fac_{max}, \max\left\{fac_{min}, fac_{scale} \cdot e ^{-1/p}\right\}\right\}.\f]
   *
   */
  template <typename RealT>
  StepControl<RealT> AdaptiveStep<RealT>::nextStep(RealT err, StepControl<RealT> prev_step, uint8_t method_order)
  {
    StepControl<RealT> next_step = prev_step;

    double h_mult = std::min(params_.fac_max_, std::max(params_.fac_scale_ * std::pow(err, GridKit::MINUS_ONE<RealT> / method_order), params_.fac_min_));

    next_step.accept_     = err <= 1;
    next_step.step_size_ *= h_mult;

    return next_step;
  }

  template class AdaptiveStep<double>;
} // namespace Integrator
