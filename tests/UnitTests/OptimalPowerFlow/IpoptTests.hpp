#pragma once

#include <limits>

#include <IpIpoptApplication.hpp>

#include <GridKit/Solver/Optimization/IpoptInterface.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "SystemTestData.hpp"
#include "TestHelpers.hpp"

namespace GridKit
{
  namespace Testing
  {
    class IpoptTests
    {
      using RealT   = double;
      using IdxT    = size_t;
      using SystemT = OPF::SystemModel<RealT, IdxT>;

    public:
      TestOutcome reference()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= optimize(system);
        if (!success)
        {
          return success.report(__func__);
        }

        const RealT* x = system.x().getData();

        success *= isEqual(x[1], ZERO<RealT>, BALANCE_TOL);
        success *= !isEqual(x[3], ZERO<RealT>, BALANCE_TOL) && !isEqual(x[5], ZERO<RealT>, BALANCE_TOL);

        return success.report(__func__);
      }

      TestOutcome solve()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= optimize(system);
        if (!success)
        {
          return success.report(__func__);
        }

        success *= system.evaluateEqualities() == 0;
        success *= system.evaluateInequalities() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* c = system.c().getData();
        const RealT* d = system.d().getData();
        for (const IdxT i : SystemTestData::BALANCE_ROWS)
        {
          success *= isEqual(c[i], ZERO<RealT>, BALANCE_TOL);
        }
        success *= isEqual(d[SystemTestData::FLOW_ROWS[0]], SystemTestData::SMAX * SystemTestData::SMAX, BALANCE_TOL);

        const RealT* x = system.x().getData();

        success *= x[6] > x[8];
        success *= x[8] > ZERO<RealT>;

        return success.report(__func__);
      }

    private:
      static constexpr RealT BALANCE_TOL = 1.0e4 * std::numeric_limits<RealT>::epsilon();

      static bool optimize(SystemT& system)
      {
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetIntegerValue("print_level", 0);
        app->Options()->SetNumericValue("bound_relax_factor", 0.0);
        app->Options()->SetNumericValue("tol", 1.0e-10);
        if (app->Initialize() != Ipopt::Solve_Succeeded)
        {
          return false;
        }

        Ipopt::SmartPtr<Ipopt::TNLP> problem = new AnalysisManager::SteadyState::IpoptInterface<RealT, IdxT>(&system);
        return app->OptimizeTNLP(problem) == Ipopt::Solve_Succeeded;
      }
    };
  } // namespace Testing
} // namespace GridKit
