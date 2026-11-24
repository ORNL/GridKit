/** @file CurrentCircle.hpp
 *  @brief Off-axis current capacity for WECC renewable electrical controls.
 */
#pragma once

#include <cmath>
#include <limits>

#include <GridKit/CommonMath.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Smooth square root of Imax^2 - high^2 with finite slope at the circle.
    template <typename RealT, typename T>
    [[gnu::always_inline]] inline T currentCircle(RealT limit, T high)
    {
      const T     x          = ((HALF<RealT> * limit - HALF<RealT> * high) * (HALF<RealT> * limit + HALF<RealT> * high)) * FOUR<RealT>;
      const RealT root_width = ONE<RealT> / Math::MU<RealT>;
      const RealT knee       = 100.0 * std::numeric_limits<RealT>::epsilon() / Math::MU<RealT>;

      const T absolute    = std::abs(x);
      const T normalizer  = absolute + knee;
      const T scaled_x    = absolute / normalizer;
      const T scaled_knee = knee / normalizer;
      const T magnitude   = normalizer * std::sqrt(scaled_x * scaled_x + scaled_knee * scaled_knee);
      const T half_sum    = HALF<RealT> * magnitude + HALF<RealT> * absolute;
      const T conjugate   = QUARTER<RealT> * knee * (knee / half_sum);
      const T hinged      = HALF<RealT> * x + HALF<RealT> * absolute + conjugate;
      // Dividing by the root keeps a large finite circle from overflowing
      return hinged / std::sqrt(hinged + root_width * root_width);
    }
  } // namespace EMT
} // namespace GridKit
