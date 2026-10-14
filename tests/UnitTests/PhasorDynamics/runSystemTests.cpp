#include "SystemTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.constructor();
  result += test.composer();
  result += test.reallocateAfterTopologyChange();
  result += test.modelVectorsAliasSystemStorage();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobian();
#endif

  result += test.allocationError();
  result += test.ieeestConstruction<0>();
  result += test.ieeestConstruction<1>();
  result += test.ieeestConstruction<2>();
  result += test.ieeestConstruction<3>();
  result += test.ieeestConstruction<4>();
  result += test.componentInitializationError();
  result += test.signalError();

  return result.summary();
}
