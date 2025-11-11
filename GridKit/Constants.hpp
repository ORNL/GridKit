#pragma once

#include <limits>

namespace GridKit
{
  template <typename IdxT>
  inline constexpr IdxT INVALID_INDEX = std::numeric_limits<IdxT>::max();

  template <typename RealT>
  inline constexpr RealT ZERO      = 0.0;

  template <typename RealT>
  inline constexpr RealT ONE       = 1.0;

  template <typename RealT>
  inline constexpr RealT TWO       = 2.0;

  template <typename RealT>
  inline constexpr RealT THREE     = 3.0;

  template <typename RealT>
  inline constexpr RealT HALF      = 0.5;

  template <typename RealT>
  inline constexpr RealT QUARTER   = 0.25;

  template <typename RealT>
  inline constexpr RealT MINUS_ONE = -1.0;
} // namespace GridKit
