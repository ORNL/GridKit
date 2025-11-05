#include <iostream>

#include <GridKit/Utilities/Testing.hpp>

/**
 * @brief Standalone example that computes the derivative of a scalar function
 * by automatic differentiation via Enzyme.
 *
 * TODO: Convert this into a unit test.
 */

double square(double x)
{
  return x * x;
}

double __enzyme_autodiff(double (*)(double), ...);

double dsquare(double x)
{
  return __enzyme_autodiff(square, x);
}

int main()
{
  int    fail = 0;
  double var  = 5.0;
  double sq   = square(var);
  double dsq  = dsquare(var);
  std::cout << "x = " << var << ", x^2 = " << sq << ", d(x^2)/dx = " << dsq << "\n";
  if (!GridKit::Testing::isEqual(dsq, 2.0 * var))
  {
    fail++;
    std::cout << "Result incorrect\n";
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
