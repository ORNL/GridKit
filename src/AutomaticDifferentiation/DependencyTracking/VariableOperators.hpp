/**
    \file VariableOperators.hpp

*/
#pragma once

namespace GridKit
{
  namespace DependencyTracking
  {
    //------------------------------------
    // non-member operators and functions
    //------------------------------------

    // unary -
    const Variable operator-(const Variable& v)
    {
      return -1.0 * v;
    }

    // +
    const Variable operator+(const Variable& lhs, const Variable& rhs)
    {
      return Variable(lhs) += rhs;
    }

    const Variable operator+(const Variable& lhs, const double& rhs)
    {
      return Variable(lhs) += rhs;
    }

    const Variable operator+(const double& lhs, const Variable& rhs)
    {
      return Variable(rhs) += lhs;
    }

    // -
    const Variable operator-(const Variable& lhs, const Variable& rhs)
    {
      return Variable(lhs) -= rhs;
    }

    const Variable operator-(const Variable& lhs, const double& rhs)
    {
      return Variable(lhs) -= rhs;
    }

    const Variable operator-(const double& lhs, const Variable& rhs)
    {
      return Variable(lhs) -= rhs;
    }

    // *
    const Variable operator*(const Variable& lhs, const Variable& rhs)
    {
      return Variable(lhs) *= rhs;
    }

    const Variable operator*(const Variable& lhs, const double& rhs)
    {
      return Variable(lhs) *= rhs;
    }

    const Variable operator*(const double& lhs, const Variable& rhs)
    {
      return Variable(lhs) *= rhs;
    }

    // /
    const Variable operator/(const Variable& lhs, const Variable& rhs)
    {
      return Variable(lhs) /= rhs;
    }

    const Variable operator/(const Variable& lhs, const double& rhs)
    {
      return Variable(lhs) /= rhs;
    }

    const Variable operator/(const double& lhs, const Variable& rhs)
    {
      return Variable(lhs) /= rhs;
    }

    // ==
    bool operator==(const Variable& lhs, const Variable& rhs)
    {
      return lhs.getValue() == rhs.getValue();
    }

    bool operator==(const Variable& lhs, const double& rhs)
    {
      return lhs.getValue() == rhs;
    }

    bool operator==(const double& lhs, const Variable& rhs)
    {
      return lhs == rhs.getValue();
    }

    // !=
    bool operator!=(const Variable& lhs, const Variable& rhs)
    {
      return !(lhs.getValue() == rhs.getValue());
    }

    bool operator!=(const Variable& lhs, const double& rhs)
    {
      return !(lhs.getValue() == rhs);
    }

    bool operator!=(const double& lhs, const Variable& rhs)
    {
      return !(lhs == rhs.getValue());
    }

    // <
    bool operator<(const Variable& lhs, const Variable& rhs)
    {
      return lhs.getValue() < rhs.getValue();
    }

    bool operator<(const Variable& lhs, const double& rhs)
    {
      return lhs.getValue() < rhs;
    }

    bool operator<(const double& lhs, const Variable& rhs)
    {
      return lhs < rhs.getValue();
    }

    // >
    bool operator>(const Variable& lhs, const Variable& rhs)
    {
      return lhs.getValue() > rhs.getValue();
    }

    bool operator>(const Variable& lhs, const double& rhs)
    {
      return lhs.getValue() > rhs;
    }

    bool operator>(const double& lhs, const Variable& rhs)
    {
      return lhs > rhs.getValue();
    }

    // <=
    bool operator<=(const Variable& lhs, const Variable& rhs)
    {
      return lhs < rhs || lhs == rhs;
    }

    bool operator<=(const Variable& lhs, const double& rhs)
    {
      return lhs < rhs || lhs == rhs;
    }

    bool operator<=(const double& lhs, const Variable& rhs)
    {
      return lhs < rhs || lhs == rhs;
    }

    // >=
    bool operator>=(const Variable& lhs, const Variable& rhs)
    {
      return lhs > rhs || lhs == rhs;
    }

    bool operator>=(const Variable& lhs, const double& rhs)
    {
      return lhs > rhs || lhs == rhs;
    }

    bool operator>=(const double& lhs, const Variable& rhs)
    {
      return lhs > rhs || lhs == rhs;
    }

    //------------------------------------
    // non-member operators
    //------------------------------------

    /**
        @brief Stream insertion operator for variables.
    */
    std::ostream& operator<<(std::ostream& os, const Variable& v)
    {
      v.print(os);
      return os;
    }

    /**
        @brief Adds all dependencies of a variable, i.e. v << u means
        that v depends on u.

        \attention Use only if you know what you are doing!
    */
    Variable& operator<<(Variable& u, const Variable& v)
    {
      u.addDependencies(v);
      return u;
    }

    /**
        @brief Stream extraction operator for variables.
    */
    std::istream& operator>>(std::istream& is, Variable& v)
    {
      return is >> v();
    }

    //-------------------------------------------
    // definitions of cmath functions derivatives
    //-------------------------------------------

    /// Derivative of sine
    inline double sin_derivative(double x)
    {
      return std::cos(x);
    }

    /// Derivative of cosine
    inline double cos_derivative(double x)
    {
      return -std::sin(x);
    }

    /// Derivative of tangent
    inline double tan_derivative(double x)
    {
      return 1.0 / (std::cos(x) * std::cos(x));
    }

    /// Derivative of arc sine
    inline double asin_derivative(double x)
    {
      return std::sqrt(1.0 - x * x);
    }

    /// Derivative of arc cosine
    inline double acos_derivative(double x)
    {
      return -std::sqrt(1.0 - x * x);
    }

    /// Derivative of arc tangent
    inline double atan_derivative(double x)
    {
      return 1.0 / (1.0 + x * x);
    }

    /// Derivative of quadrant-based arc tangent with respect to first argument
    inline double atan2_derivative1(double y, double x)
    {
      return atan_derivative(y / x) / x;
    }

    /// Derivative of quadrant-based arc tangent with respect to second argument
    inline double atan2_derivative2(double y, double x)
    {
      return -atan_derivative(y / x) * y * y / x;
    }

    /// Derivative of hyperbolic sine
    inline double sinh_derivative(double x)
    {
      return std::cosh(x);
    }

    /// Derivative of hyperbolic cosine
    inline double cosh_derivative(double x)
    {
      return std::sinh(x);
    }

    /// Derivative of hyperbolic tangent
    inline double tanh_derivative(double x)
    {
      return 1.0 / (std::cosh(x) * std::cosh(x));
    }

    /// Derivative of hyperbolic arc sine
    inline double asinh_derivative(double x)
    {
      return std::sqrt(1.0 + x * x);
    }

    /// Derivative of hyperbolic arc cosine
    inline double acosh_derivative(double x)
    {
      return std::sqrt(x * x - 1.0);
    }

    /// Derivative of hyperbolic arc tangent
    inline double atanh_derivative(double x)
    {
      return 1.0 / (1.0 - x * x);
    }

    /// Derivative of exponential function.
    inline double exp_derivative(double x)
    {
      return std::exp(x);
    }

    /// Derivative of natural logarithm function.
    inline double log_derivative(double x)
    {
      return 1.0 / x;
    }

    /// Derivative of logarithm to base 10 function: 1/(x*log(10)).
    inline double log10_derivative(double x)
    {
      return 0.4342944819032518 / x;
    }

    /// Derivative of square root function.
    inline double sqrt_derivative(double x)
    {
      return 0.5 / std::sqrt(x);
    }

    /// Derivative of absolute value function.
    inline double abs_derivative(double x)
    {
      return x > 0.0 ? 1.0 : -1.0;
    }

  } // namespace DependencyTracking
} // namespace GridKit

// Add all mathematical functions to the namespace std so that,
// for instance, std::sin(x) also works with objects of type Variable.
namespace std
{

#define IMPL_FUN_1(FUN, DER)                                                                       \
  inline GridKit::DependencyTracking::Variable FUN(const GridKit::DependencyTracking::Variable& x) \
  {                                                                                                \
    double                                val = FUN(x());                                          \
    double                                der = DER(x());                                          \
    GridKit::DependencyTracking::Variable res(x); /* copy derivatives of x*/                       \
    res.setValue(val);                            /* set function value f(x) */                    \
    res.scaleDependencies(der);                   /* compute derivatives of f(x) */                \
    return res;                                                                                    \
  }

#define IMPL_FUN_2(FUN, DER1, DER2)                                                                         \
  inline GridKit::DependencyTracking::Variable FUN(const GridKit::DependencyTracking::Variable& x,          \
                                                   const GridKit::DependencyTracking::Variable& y)          \
  {                                                                                                         \
    double                                val  = FUN(x(), y());                                             \
    double                                der1 = DER1(x(), y());                                            \
    double                                der2 = DER2(x(), y());                                            \
    GridKit::DependencyTracking::Variable res1(x);   /* copy derivatives of x*/                             \
    GridKit::DependencyTracking::Variable res2(y);   /* copy derivatives of y*/                             \
    res1.setValue(val);                              /* set function value f(x, y) */                       \
    res2.setValue(val);                              /* set function value f(x, y) */                       \
    res1.scaleDependencies(der1);                    /* compute derivatives of f(x, y) with respect to x */ \
    res2.scaleDependencies(der2);                    /* compute derivatives of f(x, y) with respect to y */ \
    GridKit::DependencyTracking::Variable res(res1); /* add derivatives with respect to x to return */      \
    res.addDependencies(res2);                       /* add derivatives with respect to y to return */      \
    return res;                                                                                             \
  }

  IMPL_FUN_1(sin, GridKit::DependencyTracking::sin_derivative)
  IMPL_FUN_1(cos, GridKit::DependencyTracking::cos_derivative)
  IMPL_FUN_1(tan, GridKit::DependencyTracking::tan_derivative)
  IMPL_FUN_1(asin, GridKit::DependencyTracking::asin_derivative)
  IMPL_FUN_1(acos, GridKit::DependencyTracking::acos_derivative)
  IMPL_FUN_1(atan, GridKit::DependencyTracking::atan_derivative)
  IMPL_FUN_1(sinh, GridKit::DependencyTracking::sinh_derivative)
  IMPL_FUN_1(cosh, GridKit::DependencyTracking::cosh_derivative)
  IMPL_FUN_1(tanh, GridKit::DependencyTracking::tanh_derivative)
  IMPL_FUN_1(exp, GridKit::DependencyTracking::exp_derivative)
  IMPL_FUN_1(log, GridKit::DependencyTracking::log_derivative)
  IMPL_FUN_1(log10, GridKit::DependencyTracking::log10_derivative)
  IMPL_FUN_1(sqrt, GridKit::DependencyTracking::sqrt_derivative)
  IMPL_FUN_1(abs, GridKit::DependencyTracking::abs_derivative)

  IMPL_FUN_2(atan2, GridKit::DependencyTracking::atan2_derivative1, GridKit::DependencyTracking::atan2_derivative2)

#undef IMPL_FUN_1
#undef IMPL_FUN_2

} // namespace std
