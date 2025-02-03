#include <iostream>
#include <limits>

double square(double x) {
  return x * x;
}

double __enzyme_autodiff(double(*)(double), ...);
double dsquare(double x) {
  return __enzyme_autodiff(square, x);
}

int main()
{
  int fail = 0;
  double var = 5.0;
  double sq  = square(var);
  double dsq = dsquare(var);
  std::cout << "x = " << var << ", x^2 = " << sq << ", d(x^2)/dx = " << dsq << "\n"; 
  if (std::abs(dsq - 2.0*var) > std::numeric_limits<double>::epsilon())
  {
    fail++;
    std::cout << "Result incorrect\n";
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
