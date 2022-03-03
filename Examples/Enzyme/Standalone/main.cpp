#include <iostream>

int square(int x) {
  return x * x;
}

int __enzyme_autodiff(int(*)(int), ...);
int dsquare(int x) {
  return __enzyme_autodiff(square, x);
}

int main() {
  int fail = 0;
  if (dsquare(5) != 25)
    fail++;
  return fail;
}
