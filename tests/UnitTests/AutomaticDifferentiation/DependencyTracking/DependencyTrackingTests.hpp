#include <iomanip>
#include <iostream>

#include <AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class SparsityPatternTests
    {
    public:
      SparsityPatternTests()  = default;
      ~SparsityPatternTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome variable()
      {
        TestStatus success = true;

        const size_t                              n = 3;
        std::vector<DependencyTracking::Variable> x(n), p(n), f(n);

        // decide x, y, and z are variables
        for (size_t i = 0; i < n; ++i)
          x[i].setVariableNumber(i);

        // ...and sigma, rho, and beta are constant parameters
        for (size_t i = 0; i < n; ++i)
          p[i].setFixed(true);

        // initialize independent variables
        x[0]() = 8.0;
        x[1]() = 20.0;
        x[2]() = 2.0 / 3.0;

        // set constant parameter values
        p[0] = 10.0;
        p[1] = 8.0 / 3.0;
        p[2] = 28.0;

        // The residualFunction computes f
        residualFunction(f, x, p);

        // Check dependenices of f[0] (depends on x[0] and x[1])
        {
          const DependencyTracking::Variable::DependencyMap& dependencies =
              (f[0]).getDependencies();

          success *= (dependencies.size() == 2);
          success *= (dependencies.find(0) != dependencies.end());
          success *= (dependencies.find(1) != dependencies.end());
          success *= (dependencies.find(2) == dependencies.end());
        }

        // Check dependencies of f[2] (depends on x[0], x[1], x[2])
        {
          const DependencyTracking::Variable::DependencyMap& dependencies =
              (f[2]).getDependencies();

          success *= (dependencies.size() == 3);
          success *= (dependencies.find(0) != dependencies.end());
          success *= (dependencies.find(1) != dependencies.end());
          success *= (dependencies.find(2) != dependencies.end());
        }

        return success.report(__func__);
      }

    private:
      template <typename T>
      void residualFunction(std::vector<T>&       f,
                            const std::vector<T>& x,
                            const std::vector<T>& p)
      {
        const T u = x[0] * x[1];

        f[0] = p[0] * (x[1] - x[0]);        // sigma*(y - x)
        f[1] = x[0] * (p[1] - x[2]) - x[1]; // x*(rho - z) - y
        f[2] = u - p[2] * x[2];             // x*y - beta*z
      }
    };

  } // namespace Testing
} // namespace GridKit
