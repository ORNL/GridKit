#pragma once

#include <cmath>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Math
  {
    /**
     * @brief Scaled sigmoid activation function
     *
     * @note The sigmoid constant (mu) value of the constant to balance accuracy
     * and finite derivatives. Large values more closely approximate a step
     * function, but lead to inf or NaN derivatives. The value of 240 corresponds
     * to a time step of 1/4 of a 60Hz cycle.
     *
     * @tparam ScalarT - scalar data type
     *
     * @param[in] x
     * @return value of the sigmoid function
     */
    template <class ScalarT>
    __attribute__((always_inline)) inline ScalarT sigmoid(const ScalarT x)
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;
      static constexpr RealT MU = 240.0;
      return ONE<RealT> / (ONE<RealT> + std::exp(-MU * x));
    }

    /**
     * @brief Low indicator function for regulator limits
     *
     * @tparam ScalarT - Scalar data type
     * @tparam RealT - Real data type (see GridKit::ScalarTraits<ScalarT>::RealT)
     *
     * @param[in] limit_min - Minimum limit
     * @param[in] x - State variable
     * @param[in] f - Conditional derivative of state variable
     * @return Scalar value indicating limit activation
     */
    template <class ScalarT, typename RealT>
    __attribute__((always_inline)) inline ScalarT indicator_low(
        const RealT limit_min, 
        const ScalarT x, 
        const ScalarT f)
    {
      return sigmoid(limit_min - x) * sigmoid(-f);
    }

    /**
     * @brief High indicator function for regulator limits
     *
     * @tparam ScalarT - Scalar data type
     * @tparam RealT - Real data type (see GridKit::ScalarTraits<ScalarT>::RealT)
     *
     * @param[in] limit_max - Maximum limit
     * @param[in] x - State variable
     * @param[in] f - Conditional derivative of state variable
     * @return Scalar value indicating limit activation
     */
    template <class ScalarT, typename RealT>
    __attribute__((always_inline)) inline ScalarT indicator_high(
        const RealT limit_max, 
        const ScalarT x, 
        const ScalarT f)
    {
      return sigmoid(x - limit_max) * sigmoid(f);
    }

    /**
     * @brief Net Indicator function for regulator limits
     *
     * @tparam ScalarT - Scalar data type
     * @tparam RealT - Real data type (see GridKit::ScalarTraits<ScalarT>::RealT)
     *
     * @param[in] limit_min - Minimum limit
     * @param[in] limit_max - Maximum limit
     * @param[in] x - State variable
     * @param[in] f - Conditional derivative of state variable
     * @return Scalar value indicating limit activation
     */
    template <class ScalarT, typename RealT>
    __attribute__((always_inline)) inline ScalarT indicator(
        const RealT limit_min, 
        const RealT limit_max, 
        const ScalarT x, 
        const ScalarT f)
    {
      return (ONE<RealT> - indicator_low(limit_min, x, f)) * (ONE<RealT> - indicator_high(limit_max, x, f));
    }
  } // namespace Math
} // namespace GridKit
