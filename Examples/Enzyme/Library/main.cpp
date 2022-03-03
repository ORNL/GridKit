#include "library.hpp"

int __enzyme_autodiff(int (*)(int), ...);

int main() {
  int fail = 0;
  if (__enzyme_autodiff(square, 5) != 25)
    fail++;
  return fail;
}
