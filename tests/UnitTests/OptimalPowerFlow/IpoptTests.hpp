#pragma once

#include <limits>
#include <vector>

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
      // Ipopt stops with the constraints and the stationarity within its tolerance. They reach 5e-13 and 2e-13,
      // and the marginal prices 8e-13.
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

        success *= isEqual(system.getBus(1)->Vi(), ZERO<RealT>, kTol);
        success *= !isEqual(system.getBus(2)->Vi(), ZERO<RealT>, kTol);
        success *= !isEqual(system.getBus(3)->Vi(), ZERO<RealT>, kTol);

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

        const RealT* g = system.g().getData();
        for (const IdxT i : SystemTestData::BALANCE_ROWS)
        {
          success *= isEqual(g[i], ZERO<RealT>, kTol);
        }
        success *= isEqual(g[SystemTestData::FLOW_ROWS[0]], SystemTestData::SMAX * SystemTestData::SMAX, kTol);

        const RealT p1 = system.getGenerator("gen_1")->P();
        const RealT p3 = system.getGenerator("gen_3")->P();

        success *= p1 > p3;
        success *= p3 > ZERO<RealT>;

        return success.report(__func__);
      }

      /// The negative of the active power multiplier of a bus is the marginal
      /// cost of a generator there within its limits
      TestOutcome marginalPrices()
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

        const auto price = [&system](IdxT number)
        {
          const auto P = static_cast<size_t>(OPF::BusInternalConstraints::P);
          return -system.getBus(number)->lambda().getData()[P];
        };
        const auto marginalCost = [&system](const char* id, RealT c1)
        { return c1 + TWO<RealT> * SystemTestData::C2 * system.getGenerator(id)->P(); };

        success *= isEqual(price(1), marginalCost("gen_1", SystemTestData::C1_GEN_1), kTol);
        success *= isEqual(price(3), marginalCost("gen_3", SystemTestData::C1_GEN_3), kTol);

        return success.report(__func__);
      }

      /// The multipliers left in the model make the Lagrangian stationary,
      /// \f$\nabla f + J^T \lambda - z_l + z_u = 0\f$
      TestOutcome multipliers()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= optimize(system);
        success *= system.evaluateGradient() == 0 && system.evaluateJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        std::vector<RealT> stationarity(system.gradient().getData(), system.gradient().getData() + system.size());
        auto*              jacobian = system.getCsrJacobian();
        const IdxT*        row_ptrs = jacobian->getRowData();
        const IdxT*        cols     = jacobian->getColData();
        const RealT*       values   = jacobian->getValues();
        const RealT*       lambda   = system.lambda().getData();
        for (IdxT row = 0; row < jacobian->getNumRows(); ++row)
        {
          for (IdxT k = row_ptrs[row]; k < row_ptrs[row + 1]; ++k)
          {
            stationarity[cols[k]] += values[k] * lambda[row];
          }
        }

        const RealT* z_lower = system.zLower().getData();
        const RealT* z_upper = system.zUpper().getData();
        for (IdxT j = 0; j < system.size(); ++j)
        {
          success *= isEqual(stationarity[j] - z_lower[j] + z_upper[j], ZERO<RealT>, kTol);
        }

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
