#include "GridKit/Solver/Dynamic/Rosenbrock.hpp"
#include "RosenbrockTests.hpp"

int main()
{
  using namespace GridKit;
  using namespace GridKit::Testing;

  GridKit::Testing::TestingResults               result;
  GridKit::Testing::RosenbrockTests<double, int> test;

  result += test.test_order(Integrator::Rosenbrock<double, int>::Tableau::lin_implicit_euler());

  return result.summary();
}
