#include "SystemTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults              result;
  GridKit::Testing::SystemTests<double, size_t> test;

  result += test.constructor();
  result += test.composer();
  result += test.residualAssemblyIsIdempotent();
  result += test.networkAdmittanceMergesDuplicateStamps();
  result += test.networkAdmittanceTracksParameterChanges();
  result += test.networkAdmittancePreservesFaultEvents();
  result += test.reallocateAfterTopologyChange();
  result += test.modelVectorsAliasSystemStorage();
  result += test.componentsShareEvaluationContext();
#ifdef GRIDKIT_ENABLE_ENZYME
  result += test.jacobianAssemblyTracksCachedContributions();
  result += test.jacobian();
#endif

  result += test.allocationError();
  result += test.componentInitializationError();
  result += test.signalError();

  return result.summary();
}
