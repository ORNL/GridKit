#include <iostream>
#include "library.hpp"

int __enzyme_autodiff(int (*)(int), ...);

int main() {
  int fail = 0;
  int sq = square(5);
  int dersq = __enzyme_autodiff(square, 5);

  std::cout << "x = 5, x^2 = " << sq << ", d(x^2)/dx = " << dersq << "\n"; 
  if (__enzyme_autodiff(square, 5) != 10)
  {
    --fail;
    // std::cout << "Incorrect result\n";
  }
  if(fail < 0)
    std::cout << "Result incorrect\n\n";
  else
    std::cout << "Bingo!!\n\n";
  std:: cout << fail << "\n";
  return fail;
}
