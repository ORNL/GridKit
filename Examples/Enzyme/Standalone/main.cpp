#include <iostream>

int foo(int x) {
  return x * x;
}

int __enzyme_autodiff(int(*)(int), ...);
int dfoo(int x) {
  return __enzyme_autodiff(foo, x);
}

int main() {
  int fail = 0;
  if (dfoo(5) != 25)
    fail++;
  return fail;
}
