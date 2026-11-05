/**
 * @file runModelTests.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Runs the EMT model tests.
 */

#include "ModelTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;
  GridKit::Testing::ModelTests     test;

  result += test.loadZ();
  result += test.conjugatePoles();
  result += test.invalidVectorFit();
  result += test.invalidInput();
  result += test.absoluteTolerances();

  return result.summary();
}
