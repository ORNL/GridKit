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
     * @note The sigmoid constant (mu) value is chosen to balance accuracy
     * and finite derivatives. Large values more closely approximate a step
     * function, but lead to inf or NaN derivatives.
     *
     * @tparam ScalarT - scalar data type
     *
     * @param[in] x - expected to be of order 1
     * @return value of the sigmoid function
     */
    template <class ScalarT>
    __attribute__((always_inline)) inline ScalarT sigmoid(const ScalarT x)
    {
      using RealT               = typename GridKit::ScalarTraits<ScalarT>::RealT;
      static constexpr RealT MU = 240.0;
      return ONE<RealT> / (ONE<RealT> + std::exp(-MU * x));
    }

    /**
     * @brief Derivative of the scaled sigmoid activation function
     *        (i.e., approximation to the delta dirac function)
     *
     * @tparam ScalarT - scalar data type
     *
     * @param[in] x - expected to be of order 1
     * @return value of the sigmoid function
     */
    template <class ScalarT>
    __attribute__((always_inline)) inline ScalarT dsigmoid(const ScalarT x)
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;
      return FOUR<RealT> * sigmoid(x) * (ONE<RealT> - sigmoid(x));
    }

    /**
     * @brief Low indicator function for regulator limits
     *
     * @tparam ScalarT - Scalar data type
     * @tparam RealT - Real data type (see GridKit::ScalarTraits<ScalarT>::RealT)
     *
     * @param[in] limit_min - Minimum limit
     * @param[in] x - State variable
     * @return Scalar value indicating limit activation
     */
    template <class ScalarT, typename RealT>
    __attribute__((always_inline)) inline ScalarT indicator_low(
        const RealT   limit_min,
        const ScalarT x)
    {
      return sigmoid(x - limit_min);
    }

    /**
     * @brief High indicator function for regulator limits
     *
     * @tparam ScalarT - Scalar data type
     * @tparam RealT - Real data type (see GridKit::ScalarTraits<ScalarT>::RealT)
     *
     * @param[in] limit_max - Maximum limit
     * @param[in] x - State variable
     * @return Scalar value indicating limit activation
     */
    template <class ScalarT, typename RealT>
    __attribute__((always_inline)) inline ScalarT indicator_high(
        const RealT   limit_max,
        const ScalarT x)
    {
      return sigmoid(limit_max - x);
    }

    /**
     * @brief Zero indicator function for regulator limits
     *
     * @tparam ScalarT - Scalar data type
     * @tparam RealT - Real data type (see GridKit::ScalarTraits<ScalarT>::RealT)
     *
     * @param[in] limit_max - Maximum limit
     * @param[in] x - State variable
     * @return Scalar value indicating limit activation
     */
    template <class ScalarT, typename RealT>
    __attribute__((always_inline)) inline ScalarT indicator_zero(
        const RealT   limit_min,
        const RealT   limit_max,
        const ScalarT x)
    {
      return indicator_low(limit_min, x) + indicator_high(limit_max, x) - ONE<RealT>;
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
        const RealT   limit_min,
        const RealT   limit_max,
        const ScalarT x,
        const ScalarT f)
    {
      return indicator_zero(limit_min, limit_max, f) * dsigmoid(f) + //
             (indicator_low(limit_min, x) * sigmoid(-f) +            //
              indicator_high(limit_max, x) * sigmoid(f))
                 * (ONE<RealT> - dsigmoid(f));
    }
  } // namespace Math
} // namespace GridKit
