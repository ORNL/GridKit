#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>

#include "RosenbrockTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;
  using namespace AnalysisManager::NativeDynamicSolver;

  GridKit::Testing::TestingResults               result;
  GridKit::Testing::RosenbrockTests<double, int> test;

  result += test.test_order(Rosenbrock<double, int>::Tableau::linImplicitEuler(), -1.0, -5.0);
  result += test.test_order(Rosenbrock<double, int>::Tableau::rodas5p(), -0.5, -2.5);

  return result.summary();
}
