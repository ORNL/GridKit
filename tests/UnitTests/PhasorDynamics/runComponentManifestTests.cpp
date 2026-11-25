#include "ComponentManifestTests.hpp"

int main()
{
  GridKit::Testing::TestingResults result;

  GridKit::Testing::ComponentManifestTests test;

  result += test.manifest();
  result += test.library();

  return result.summary();
}
