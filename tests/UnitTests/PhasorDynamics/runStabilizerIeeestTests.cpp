#include "StabilizerIeeestTests.hpp"

template <size_t order>
void runOrder(GridKit::Testing::TestingResults& result)
{
  std::cout << "IEEEST order " << order << '\n';
  GridKit::Testing::StabilizerIeeestTests<double, size_t, order> test;
  result += test.constructor();
  result += test.validation();
  result += test.factory();
  result += test.zeroInitialResidual();
  result += test.initialization();
  result += test.transferResponse();
  result += test.limiter();
  if constexpr (order > 0)
  {
    result += test.coefficientScaling();
  }
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.lagScaling();
  result += test.jacobian();
#endif
}

int main()
{
  GridKit::Testing::TestingResults result;
  runOrder<0>(result);
  runOrder<1>(result);
  runOrder<2>(result);
  runOrder<3>(result);
  runOrder<4>(result);
  return result.summary();
}
