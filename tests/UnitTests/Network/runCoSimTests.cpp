#include "CoSimTests.hpp"

int main()
{
  using namespace GridKit::Testing;

  TestingResults                  result;
  CoSimTests<double, std::size_t> test;

  result += test.sessions();
  result += test.rejectedReplies();
  result += test.malformedMessages();
  result += test.boundedWaits();

  return result.summary();
}
