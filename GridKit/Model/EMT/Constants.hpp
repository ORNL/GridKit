/**
 * @file Constants.hpp
 *
 * @brief Shared physical constants for EMT models.
 *
 */

#pragma once

#include <numbers>

namespace GridKit
{
  namespace EMT
  {
    namespace Constants
    {
      template <typename ScalarT>
      constexpr ScalarT pi()
      {
        return static_cast<ScalarT>(std::numbers::pi_v<double>);
      }

      template <typename ScalarT>
      constexpr ScalarT mu0()
      {
        return static_cast<ScalarT>(4.0e-7) * pi<ScalarT>();
      }

      template <typename ScalarT>
      constexpr ScalarT epsilon0()
      {
        return static_cast<ScalarT>(8.8541878128e-12);
      }
    } // namespace Constants
  } // namespace EMT
} // namespace GridKit
