#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace GridKit::PhasorDynamics
{
  /// Positive difference with derivatives on both arguments for sparse AD.
  template <class T, class U>
  [[gnu::always_inline]] inline auto positiveDifference(T a, U b)
  {
    const auto x = a - b;
    return 0.5 * (x + std::abs(x));
  }

  /// Invert a monotone characteristic on its operating range.
  template <class F, class T>
  T invertControl(F function, T target, T low, T high)
  {
    if (!std::isfinite(target) || target < function(low) - 1e-12 || target > function(high) + 1e-12)
      return std::numeric_limits<T>::quiet_NaN();
    const T guess = std::clamp(target, low, high);
    if (std::abs(function(guess) - target) <= 1e-12)
      return guess;
    for (unsigned i = 0; i < 80; ++i)
    {
      const T mid = std::midpoint(low, high);
      if (function(mid) < target)
        low = mid;
      else
        high = mid;
    }
    return std::midpoint(low, high);
  }
} // namespace GridKit::PhasorDynamics
