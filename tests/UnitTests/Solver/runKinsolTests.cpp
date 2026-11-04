#include "KinsolTests.hpp"

int main()
{
  using namespace GridKit::Testing;

  TestingResults                   result;
  KinsolTests<double, std::size_t> test;

  result += test.nonlinearResidual();
  result += test.exceptionAndRestart();

  return result.summary();
}
