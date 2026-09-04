#include <GridKit/Definitions.hpp>

#include "ComponentConnectionTests.hpp"

int main()
{
  GridKit::Testing::TestingResults                           result;
  GridKit::Testing::ComponentConnectionTests<double, size_t> test;
  GridKit::Testing::GastPtiConnectionTests<double, size_t>   gastpti;

  result += test.genrouEsdc1a();
  result += test.genrouHygov();
  result += test.regcaRepca();
  result += test.ieeestIeeet1();
  result += gastpti.genrouGastPti();
  result += gastpti.gensalGastPti();
  result += test.regcaReecb();

  return result.summary();
}
