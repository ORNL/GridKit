#include "StabilizerIeeestTests.hpp"

template <size_t order>
void runNotchOrder(GridKit::Testing::TestingResults& result)
{
  std::cout << "IEEEST notch order " << order << '\n';

  GridKit::Testing::StabilizerIeeestTests<double, size_t, order> test;

  result += test.validation();
  result += test.factory();
  result += test.steadyState();
  result += test.transferFunction();
  if constexpr (order >= 1 && order <= 3)
  {
    result += test.factorOrdering();
  }
  result += test.limiter();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif
}

int main()
{
  GridKit::Testing::TestingResults result;

  runNotchOrder<0>(result);
  runNotchOrder<1>(result);
  runNotchOrder<2>(result);
  runNotchOrder<3>(result);
  runNotchOrder<4>(result);

  return result.summary();
}
