#pragma once

#include <limits>

namespace GridKit
{
  template <typename index_type>
  inline constexpr index_type INVALID_INDEX = std::numeric_limits<index_type>::max();

  template <typename real_type>
  inline constexpr real_type ZERO = 0.0;

  template <typename real_type>
  inline constexpr real_type ONE = 1.0;

  template <typename real_type>
  inline constexpr real_type TWO = 2.0;

  template <typename real_type>
  inline constexpr real_type THREE = 3.0;

  template <typename real_type>
  inline constexpr real_type FOUR = 4.0;

  template <typename real_type>
  inline constexpr real_type HALF = 0.5;

  template <typename real_type>
  inline constexpr real_type QUARTER = 0.25;

  template <typename real_type>
  inline constexpr real_type MINUS_ONE = -1.0;
} // namespace GridKit
