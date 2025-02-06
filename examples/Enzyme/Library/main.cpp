#include <iostream>
#include <limits>
#include "Library.hpp"

int main() {
  int fail = 0;
  Model model;
  double var = 5.0;
  model.setVariable(var);
  model.evalFunction();
  model.evalDerivative();
  double sq = model.getFunctionValue();
  double dsq = model.getDerivativeValue();

  std::cout << "x = " << var << ", x^2 = " << sq << ", d(x^2)/dx = " << dsq << "\n"; 
  if (std::abs(dsq - 2.0*var) > std::numeric_limits<double>::epsilon())
  {
    fail++;
    std::cout << "Result incorrect\n";
  }
  std::cout << "Status: " << fail << "\n";
  return fail;
}
