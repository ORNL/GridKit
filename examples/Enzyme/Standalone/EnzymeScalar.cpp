#include <iostream>
#include <limits>

/**
 * @brief Standalone example that computes the derivative of a scalar function
 * by automatic differentiation via Enzyme.
 *
 * TODO: Convert this into a unit test.
 */

template <class ScalarT>
ScalarT square(ScalarT x)
{
  return x * x;
}

template <class ScalarT>
ScalarT __enzyme_autodiff(ScalarT (*)(ScalarT), ...);

template <class ScalarT>
ScalarT dsquare(ScalarT x)
{
  return __enzyme_autodiff(square<ScalarT>, x);
}

int main()
{
  int    fail = 0;
  double var  = 5.0;
  double sq   = square<double>(var);
  double dsq  = dsquare<double>(var);
  std::cout << "x = " << var << ", x^2 = " << sq << ", d(x^2)/dx = " << dsq << "\n";
  if (std::abs(dsq - 2.0 * var) > std::numeric_limits<double>::epsilon())
  {
    fail++;
    std::cout << "Result incorrect\n";
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
