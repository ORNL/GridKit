#include "IeeestTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::IeeestTests    test;
  result += test.phasorParityAndInitialization();
  result += test.transferFunction();
  result += test.jacobianAndCutout();
  result += test.validation();
  return result.summary();
}
