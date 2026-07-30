#pragma once

#include <cassert>
#include <cmath>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Math
  {
    /**
     * @brief Smoothing scale shared by CommonMath primitives
     *
     * Used by @ref sigmoid, @ref ramp, and functions composed from them to set
     * the width of smooth transitions.
     *
     * @tparam real_type - real data type
     */
    template <typename real_type>
    inline constexpr real_type MU = 240.0;

    /**
     * @brief Scaled sigmoid activation function
     *
     * @note The sigmoid constant (mu) value is chosen to balance accuracy
     * and finite derivatives. Large values more closely approximate a step
     * function, but can make the transition numerically stiff.
     *
     * @tparam scalar_type - scalar data type
     *
     * @param[in] x - expected to be of order 1
     * @return value of the sigmoid function
     */
    template <typename scalar_type>
    __attribute__((always_inline)) inline scalar_type sigmoid(const scalar_type x)
    {
      using ScalarT = scalar_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;
      return HALF<RealT> * (ONE<RealT> + std::tanh(HALF<RealT> * MU<RealT> * x));
    }

    /**
     * @brief Smooth one-sided ramp function
     *
     * Smooth approximation to max(x, 0), using a stable softplus form with
     * the same scale as the rest of CommonMath.
     *
     * @tparam scalar_type - scalar data type
     *
     * @param[in] x - expected to be of order 1
     * @return value of the smooth ramp function
     */
    template <typename scalar_type>
    __attribute__((always_inline)) inline scalar_type ramp(const scalar_type x)
    {
      using ScalarT = scalar_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      RealT   mu = MU<RealT>;
      ScalarT a  = std::abs(mu * x);
      return HALF<RealT> * (x + a / mu) + std::log1p(std::exp(-a)) / mu;
    }

    /**
     * @brief Smooth one-sided quadratic ramp
     *
     * Smooth approximation to max(x, 0)^2 via a sigmoid-gated quadratic.
     * Used for IEEE-style quadratic saturation curves.
     *
     * @note Eventually a enzyme specialization for an exact implementation
     *       would be nice, since the piecewise definition is C^1 continuous
     *
     * @tparam scalar_type - scalar data type
     *
     * @param[in] x - input signal
     * @return value of the quadratic ramp
     */
    template <typename scalar_type>
    __attribute__((always_inline)) inline scalar_type qramp(const scalar_type x)
    {
      return x * x * sigmoid(x);
    }

    /**
     * @brief Smooth binary maximum function
     *
     * Smooth approximation to max(x, y), composed from the smooth ramp
     * function.
     *
     * @tparam left_type - scalar type of x
     * @tparam right_type - scalar type of y
     *
     * @param[in] x - First input signal
     * @param[in] y - Second input signal
     * @return Smooth maximum of x and y
     *
     * @note The two input types intentionally may differ. Model equations
     * often compare a differentiable state or signal with a plain real
     * parameter, limit, or literal bound. Keeping both template parameters
     * lets the expression promote to the differentiable scalar type without
     * forcing callers to cast every parameter.
     */
    template <typename left_type, typename right_type>
    __attribute__((always_inline)) inline auto max(
        const left_type  x,
        const right_type y)
    {
      return y + ramp(x - y);
    }

    /**
     * @brief Smooth binary minimum function
     *
     * Smooth approximation to min(x, y), composed from the smooth ramp
     * function.
     *
     * @tparam left_type - scalar type of x
     * @tparam right_type - scalar type of y
     *
     * @param[in] x - First input signal
     * @param[in] y - Second input signal
     * @return Smooth minimum of x and y
     *
     * @note The two input types intentionally may differ. Model equations
     * often compare a differentiable state or signal with a plain real
     * parameter, limit, or literal bound. Keeping both template parameters
     * lets the expression promote to the differentiable scalar type without
     * forcing callers to cast every parameter.
     */
    template <typename left_type, typename right_type>
    __attribute__((always_inline)) inline auto min(
        const left_type  x,
        const right_type y)
    {
      return x - ramp(x - y);
    }

    /**
     * @brief Smooth clamp function
     *
     * Smooth approximation to min(max(x, lower), upper), composed from the
     * smooth ramp function. Lower and upper bounds may be independent types
     * (e.g. constant Real bounds or algebraic-variable bounds).
     *
     * @tparam scalar_type - scalar data type of the input signal
     * @tparam lower_type - data type of the lower bound
     * @tparam upper_type - data type of the upper bound
     *
     * @param[in] x - expected to be of order 1
     * @param[in] lower - Lower limit
     * @param[in] upper - Upper limit
     * @return value of the smooth clamp function
     */
    template <typename scalar_type, typename lower_type, typename upper_type>
    __attribute__((always_inline)) inline auto clamp(
        const scalar_type x,
        const lower_type  lower,
        const upper_type  upper)
    {
      assert(lower <= upper);
      return lower + ramp(x - lower) - ramp(x - upper);
    }

    /**
     * @brief Smooth two-sided deadband function
     *
     * Smooth approximation to x - min(max(x, lower), upper), composed from the
     * smooth ramp function.
     *
     * @tparam scalar_type - scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - Input signal
     * @param[in] lower - Lower breakpoint
     * @param[in] upper - Upper breakpoint
     * @return Smooth deadbanded value
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type deadband(
        const scalar_type x,
        const real_type   lower,
        const real_type   upper)
    {
      assert(lower <= upper);
      return ramp(x - upper) - ramp(-(x - lower));
    }

    /**
     * @brief Smooth slew-rate limiter
     *
     * Smooth approximation to min(max(f, -rate), rate).
     *
     * @tparam scalar_type - scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] f - Pre-limit derivative or rate signal
     * @param[in] rate - Symmetric positive rate limit
     * @return Slew-rate-limited value of f
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type slew(
        const scalar_type f,
        const real_type   rate)
    {
      assert(rate >= ZERO<real_type>);
      return clamp(f, -rate, rate);
    }

    /**
     * @brief Smooth linear segment contribution
     *
     * Smooth approximation to a linear segment contribution that is zero below
     * lower, linear over [lower, upper], and saturated at height above upper.
     * Callers should supply lower < upper; height may be positive or negative.
     *
     * @tparam scalar_type - scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - Input signal
     * @param[in] lower - Lower breakpoint
     * @param[in] upper - Upper breakpoint
     * @param[in] height - Saturated value above the upper breakpoint
     * @return Smooth linear segment contribution
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type linseg(
        const scalar_type x,
        const real_type   lower,
        const real_type   upper,
        const real_type   height)
    {
      assert(lower < upper);
      return height / (upper - lower) * (ramp(x - lower) - ramp(x - upper));
    }

    /**
     * @brief Smooth above-limit indicator
     *
     * @tparam scalar_type - Scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - State variable
     * @param[in] limit_min - Minimum limit
     * @return Smooth indicator that x is above limit_min
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type above(
        const scalar_type x,
        const real_type   limit_min)
    {
      return sigmoid(x - limit_min);
    }

    /**
     * @brief Smooth below-limit indicator
     *
     * @tparam scalar_type - Scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - State variable
     * @param[in] limit_max - Maximum limit
     * @return Smooth indicator that x is below limit_max
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type below(
        const scalar_type x,
        const real_type   limit_max)
    {
      return sigmoid(limit_max - x);
    }

    /**
     * @brief Smooth inside-limits indicator
     *
     * @tparam scalar_type - Scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - State variable
     * @param[in] limit_min - Minimum limit
     * @param[in] limit_max - Maximum limit
     * @return Smooth indicator that x is inside [limit_min, limit_max]
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type inside(
        const scalar_type x,
        const real_type   limit_min,
        const real_type   limit_max)
    {
      assert(limit_min <= limit_max);
      return above(x, limit_min) + below(x, limit_max) - ONE<real_type>;
    }

    /**
     * @brief Smooth outside-limits indicator
     *
     * @tparam scalar_type - Scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - State variable
     * @param[in] limit_min - Minimum limit
     * @param[in] limit_max - Maximum limit
     * @return Smooth indicator that x is outside [limit_min, limit_max]
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type outside(
        const scalar_type x,
        const real_type   limit_min,
        const real_type   limit_max)
    {
      assert(limit_min <= limit_max);
      return below(x, limit_min) + above(x, limit_max);
    }

    /**
     * @brief Smooth anti-windup indicator for a limited state variable
     *
     * @tparam scalar_type - Scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - State variable
     * @param[in] f - Pre-limit derivative of the state variable
     * @param[in] limit_min - Minimum limit
     * @param[in] limit_max - Maximum limit
     * @return Scalar value in [0, 1]: 1 when dynamics should pass through,
     *         0 when integration should be blocked.
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type indicator(
        const scalar_type x,
        const scalar_type f,
        const real_type   limit_min,
        const real_type   limit_max)
    {
      assert(limit_min <= limit_max);
      using ScalarT = scalar_type;
      using RealT   = real_type;

      ScalarT above_min = above(x, limit_min);
      ScalarT below_max = below(x, limit_max);

      return above_min * below_max +                  //
             (ONE<RealT> - below_max) * sigmoid(-f) + //
             (ONE<RealT> - above_min) * sigmoid(f);
    }

    /**
     * @brief Smooth anti-windup limited derivative
     *
     * Applies the smooth anti-windup indicator gate to a pre-limit derivative.
     * The returned value approximates the conditional-integration rule that
     * passes interior dynamics, passes restoring motion from saturated limits,
     * and blocks motion that would push further into saturation.
     *
     * @tparam scalar_type - Scalar data type
     * @tparam real_type - Real data type (see GridKit::ScalarTraits<scalar_type>::RealT)
     *
     * @param[in] x - Limited state or limited output signal
     * @param[in] f - Pre-limit derivative
     * @param[in] limit_min - Minimum limit
     * @param[in] limit_max - Maximum limit
     * @return Smooth anti-windup limited derivative
     */
    template <typename scalar_type, typename real_type>
    __attribute__((always_inline)) inline scalar_type antiwindup(
        const scalar_type x,
        const scalar_type f,
        const real_type   limit_min,
        const real_type   limit_max)
    {
      return indicator(x, f, limit_min, limit_max) * f;
    }
  } // namespace Math
} // namespace GridKit
