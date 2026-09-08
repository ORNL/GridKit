#include "InnerCurrentControlTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                           result;
  GridKit::Testing::InnerCurrentControlTests<double, size_t> test;
  result += test.innerControl();
  return result.summary();
}
