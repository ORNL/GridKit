#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class BusFaultTests
    {
    private:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

    public:
      BusFaultTests()  = default;
      ~BusFaultTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* fault =
            new PhasorDynamics::BusFault<ScalarT, IdxT>(bus);

        success *= (fault != nullptr);

        if (fault)
        {
          delete fault;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * Verifies the residual evaluates to zero for the initial conditions
       */
      TestOutcome zeroInitialResidual(bool status = false)
      {
        TestStatus success = true;

        ScalarT Vr1{1.0}; ///< Bus real voltage
        ScalarT Vi1{1.0}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>      bus(Vr1, Vi1);
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, 0.0, 1e-3, status);
        bus.allocate();
        bus.initialize();
        fault.allocate();
        fault.initialize();
        fault.evaluateResidual();
        auto&       res      = fault.getResidual();
        const auto* res_data = res.getData();
        const auto* yp       = fault.yp().getData();

        for (size_t i = 0; i < res.getSize(); ++i)
        {
          if (!isEqual(res_data[i], 0.0))
          {
            std::cout << "Incorrect result: "
                      << yp[i] << " != 0\n";
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * A test case to verify Jacobian values
       */
      TestOutcome jacobian(bool status = false)
      {
        TestStatus success = true;

        RealT R = 0.0;
        RealT X = 1e-3;

        // Jacobian via DependencyTracking
        auto dependency_tracking_jacobian = DependencyTrackingJacobian(R, X, status);

        // Jacobian via Enzyme
        auto enzyme_jacobian = EnzymeJacobian(R, X, status);

        success *= dependency_tracking_jacobian.size() == enzyme_jacobian.size();

        const auto remove_zeros = [](auto& jacobian)
        {
          for (auto& row : jacobian)
          {
            for (auto entry = row.begin(); entry != row.end();)
            {
              if (entry->second == 0.0)
              {
                entry = row.erase(entry);
              }
              else
              {
                ++entry;
              }
            }
          }
        };
        remove_zeros(dependency_tracking_jacobian);
        remove_zeros(enzyme_jacobian);

        /// Compare DependencyTracking dependencies to Enzyme's
        const size_t rows = std::min(dependency_tracking_jacobian.size(), enzyme_jacobian.size());
        for (size_t i = 0; i < rows; ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          const RealT R, const RealT X, const bool status)
      {
        DependencyTracking::Variable Vr1{1.0}; ///< Bus-1 real voltage
        DependencyTracking::Variable Vi1{1.0}; ///< Bus-1 imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>      bus(Vr1, Vi1);
        PhasorDynamics::BusFault<DependencyTracking::Variable, IdxT> fault(&bus, R, X, status);

        bus.allocate();
        fault.allocate();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + fault.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + fault.size()); // Reset bus residual indices
        }

        bus.initialize();
        fault.initialize();

        fault.updateTime(0.0, 1.0);

        bus.evaluateResidual();
        fault.evaluateResidual();

        fault.evaluateJacobian();
        auto* model_jacobian = fault.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: BusFault DependencyTracking Jacobian\n";
        model_jacobian->print();

        auto dependencies = GridKit::Testing::MapFromCsr(model_jacobian);
        const auto& bus_residual = bus.getResidual();
        const auto internal_rows = dependencies.size();
        dependencies.resize(internal_rows + bus_residual.getSize());
        for (IdxT row = 0; row < bus_residual.getSize(); ++row)
        {
          // Merge even y and odd yp indices at alpha = 1, including bus rows.
          for (const auto& [column, value] : bus_residual.getData()[row].getDependencies())
          {
            dependencies[internal_rows + row][column / 2] += value;
          }
        }
        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          const RealT R, const RealT X, const bool status)
      {
        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>      bus(Vr1, Vi1);
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, R, X, status);

        bus.allocate();
        fault.allocate();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + fault.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + fault.size()); // Reset bus residual indices
        }

        bus.initialize();
        fault.initialize();

        fault.updateTime(0.0, 1.0);

        bus.evaluateResidual();
        fault.evaluateResidual();

        bus.evaluateJacobian();
        fault.evaluateJacobian();
        fault.constructCsr();
        auto* model_jacobian = fault.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: BusFault Enzyme Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif

    }; // class BusFaultTests

  } // namespace Testing
} // namespace GridKit
