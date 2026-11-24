#pragma once

#include <cassert>

#include <GridKit/CommonMath.hpp>

namespace GridKit::Model
{
  /**
   * @brief Smooth approximation of the REGCA `rrpwr` rate limiter.
   *
   * Limits motion that increases \f$|x|\f$ while smoothly releasing
   * restoring motion toward zero. At \f$x = 0\f$, the result is the
   * symmetric smooth slew limit applied to @p f.
   *
   * @param[in] x State whose sign determines the limited direction.
   * @param[in] f Unconstrained derivative.
   * @param[in] rate Nonnegative symmetric rate limit.
   * @return Rate-limited derivative.
   */
  template <class ScalarT, class RealT>
  [[gnu::always_inline]] inline ScalarT rrpwr(
      const ScalarT x,
      const ScalarT f,
      const RealT   rate)
  {
    assert(rate >= ZERO<RealT>);

    const ScalarT t     = std::tanh(HALF<RealT> * Math::MU<RealT> * x);
    const ScalarT abs_t = std::abs(t);
    const ScalarT w_pos = HALF<RealT> * t * (t + abs_t);
    const ScalarT w_neg = HALF<RealT> * t * (t - abs_t);

    return f + (ONE<RealT> - w_pos) * Math::ramp(-f - rate)
           - (ONE<RealT> - w_neg) * Math::ramp(f - rate);
  }

  /**
   * @brief Smooth anti-windup derivative under a moving upper bound.
   *
   * Below the bound the unconstrained derivative passes. Pinned at or
   * above the bound, the state tracks min(f, fmax), so a falling bound
   * drags the state down with it. Equivalent to fixed-bound anti-windup
   * on the gap x - xmax held below zero.
   *
   * @param[in] x State limited from above.
   * @param[in] f Unconstrained derivative of x.
   * @param[in] xmax Moving upper bound on x.
   * @param[in] fmax Derivative of the moving upper bound.
   * @return Anti-windup-limited derivative.
   */
  template <class ScalarT>
  [[gnu::always_inline]] inline ScalarT awmax(
      const ScalarT x,
      const ScalarT f,
      const ScalarT xmax,
      const ScalarT fmax)
  {
    using RealT            = typename ScalarTraits<ScalarT>::RealT;
    const ScalarT gap_rate = f - fmax;
    const ScalarT below    = Math::sigmoid(xmax - x);

    return fmax
           + (below + (ONE<RealT> - below) * Math::sigmoid(-gap_rate))
                 * gap_rate;
  }
} // namespace GridKit::Model
