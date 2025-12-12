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
    __attribute__((always_inline)) inline ScalarT sigmoid(ScalarT x)
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;
      static constexpr RealT MU = 240.0;
      return ONE<RealT> / (ONE<RealT> + std::exp(-MU * x));
    }
  } // namespace Math
} // namespace GridKit
