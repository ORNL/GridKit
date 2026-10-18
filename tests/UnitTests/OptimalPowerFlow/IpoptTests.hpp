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
      // Ipopt stops with every constraint within its tolerance. The balance rows reach 4e-13.
      static constexpr RealT kTol = static_cast<RealT>(1.0e-10);

      TestOutcome reference()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
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

        success *= isEqual(x[1], ZERO<RealT>, kTol);
        success *= !isEqual(x[3], ZERO<RealT>, kTol) && !isEqual(x[5], ZERO<RealT>, kTol);

        return success.report(__func__);
      }

      TestOutcome solve()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
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
          success *= isEqual(c[i], ZERO<RealT>, kTol);
        }
        success *= isEqual(d[SystemTestData::FLOW_ROWS[0]], SystemTestData::SMAX * SystemTestData::SMAX, kTol);

        const RealT* x = system.x().getData();

        success *= x[6] > x[8];
        success *= x[8] > ZERO<RealT>;

        return success.report(__func__);
      }

    private:
      static bool optimize(SystemT& system)
      {
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetIntegerValue("print_level", 0);
        app->Options()->SetNumericValue("bound_relax_factor", 0.0);
        app->Options()->SetNumericValue("tol", kTol);
        app->Options()->SetNumericValue("constr_viol_tol", kTol);
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
