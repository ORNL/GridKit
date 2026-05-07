#include "ExciterIeeet1Tests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ExciterIeeet1Tests<double, size_t> test;

  result += test.constructor();

  return result.summary();
}
