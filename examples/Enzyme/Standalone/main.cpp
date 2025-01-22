#include <iostream>
#include <limits>

double square(double x) {
  return x * x;
}

double __enzyme_autodiff(double(*)(double), ...);
double dsquare(double x) {
  return __enzyme_autodiff(square, x);
}

double sq  = square(5.0);
double dsq = dsquare(5.0);

int main()
{
  int fail = 0;
  std::cout << "x = 5, x^2 = " << sq << ", d(x^2)/dx = " << dsq << "\n"; 
  if (std::abs(dsq - 10.0) > std::numeric_limits<double>::epsilon())
  {
    fail++;
    std::cout << "Result incorrect\n";
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
