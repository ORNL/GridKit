#include "IdaTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults           result;
  GridKit::Testing::IdaTests<double, size_t> test;

  result += test.test();

  return result.summary();
}
