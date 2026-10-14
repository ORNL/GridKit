#include "ImplicitTrapezoidalTests.hpp"

int main()
{
  GridKit::Testing::TestingResults           result;
  GridKit::Testing::ImplicitTrapezoidalTests tests;
  result += tests.differentialAlgebraicSystem();
  return result.summary();
}
