#include <iostream>

int square(int x) {
  return x * x;
}

int __enzyme_autodiff(int(*)(int), ...);
int dsquare(int x) {
  return __enzyme_autodiff(square, x);
}

int sq  = square(5);
int dsq = dsquare(5);

int main()
{
  int fail = 0;
  std::cout << "x = 5, x^2 = " << sq << ", d(x^2)/dx = " << __enzyme_autodiff(square, 5) << "\n"; 
  if (dsq != 10)
    fail++;
  return fail;
}
