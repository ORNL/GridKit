#include <array>
#include <cmath>

#include <GridKit/Solver/SteadyState/Kinsol.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class KinsolTests
    {
      using Solver = AnalysisManager::Sundials::Kinsol<scalar_type, index_type>;

    public:
      TestOutcome nonlinearResidual()
      {
        TestStatus outcome;

        const auto residual = [](auto x, auto f)
        {
          f[0] = x[0] * x[0] + x[1] - 5.0;
          f[1] = x[0] + x[1] * x[1] - 5.0;

          return 0;
        };

        Solver solver(2);

        solver.configureSimulation();
        solver.setTolerance(1e-10);
        solver.setFunctionError(1e-8);

        std::array values{2.0, 2.5};

        solver.runSimulation(values, residual);

        const double root = (std::sqrt(21.0) - 1.0) / 2.0;

        outcome *= isEqual(values[0], root, 1e-9);
        outcome *= isEqual(values[1], root, 1e-9);

        return outcome.report(__func__);
      }

      TestOutcome exceptionAndRestart()
      {
        TestStatus outcome;

        bool       fail     = true;
        const auto residual = [&](auto x, auto f)
        {
          if (fail)
          {
            throw 42;
          }

          f[0] = x[0] - 2.0;
          f[1] = x[1] - 3.0;

          return 0;
        };

        Solver solver(2);

        solver.configureSimulation();
        solver.setTolerance(1e-10);

        std::array values{1.0, 1.0};

        outcome *= throws<int>([&]
                               { solver.runSimulation(values, residual); });

        fail = false;
        solver.runSimulation(values, residual, false);
        outcome *= isEqual(values[0], 2.0);
        outcome *= isEqual(values[1], 3.0);

        values = {4.0, 6.0};
        solver.runSimulation(values, residual);
        outcome *= isEqual(values[0], 2.0);
        outcome *= isEqual(values[1], 3.0);

        solver.runSimulation(
          values,
          [](auto x, auto f)
          {
            f[0] = x[0] - 7.0;
            f[1] = x[1] + 1.0;
            return 0;
          }
        );
        outcome *= isEqual(values[0], 7.0);
        outcome *= isEqual(values[1], -1.0);

        return outcome.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
