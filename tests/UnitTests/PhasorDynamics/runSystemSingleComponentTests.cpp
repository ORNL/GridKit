#include "SystemSingleComponentTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults                             result;
  GridKit::Testing::SystemSingleComponentTests<double, size_t> test;

  result += test.branch();
  result += test.bus();
  result += test.busFault();
  result += test.ieeet1();
  result += test.load();
  result += test.loadZIP();
  result += test.genrou();
  result += test.genClassical();
  result += test.tgov1();

  // The following components are not tested here because they require signal attachments
  // PhasorDynamics::Exciter::SexsPti
  // PhasorDynamics::Stabilizer::Ieeest

  return result.summary();
}
