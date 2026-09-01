#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
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

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]));
        }

        return success.report(__func__);
      }

      /**
       * A test case to verify a cleared fault drops out of the bus equations.
       *
       * After the fault clears it no longer injects current, so the Jacobian
       * of the bus equations must not depend on the fault variables. A stale
       * coupling left behind stalls the integrator at fault clearing.
       */
      TestOutcome jacobianAfterClearing()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>      bus(1.0, 0.0);
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, 0.0, 0.01, 0);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        system.addBus(&bus);
        system.addFault(&fault);
        system.allocate();
        system.initialize();

        // Integrator callbacks before the fault
        system.updateTime(0.0, 0.0);
        system.evaluateResidual();
        system.evaluateJacobian();

        fault.setStatus(true);

        // Integrator callbacks while the fault is on
        system.updateTime(1.0, 0.0);
        system.evaluateResidual();
        system.evaluateJacobian();

        fault.setStatus(false);

        // Integrator callbacks after clearing
        system.updateTime(2.0, 0.0);
        system.evaluateResidual();
        system.evaluateJacobian();

        std::vector<DependencyTracking::Variable::DependencyMap> jacobian =
            GridKit::Testing::MapFromCsr(system.getCsrJacobian());

        // After clearing, the bus current balance should not depend on the
        // fault current
        success *= jacobianEntryIsZero(jacobian, 0, 2); // d(bus Ir eq)/d(fault Ir)
        success *= jacobianEntryIsZero(jacobian, 0, 3); // d(bus Ir eq)/d(fault Ii)
        success *= jacobianEntryIsZero(jacobian, 1, 2); // d(bus Ii eq)/d(fault Ir)
        success *= jacobianEntryIsZero(jacobian, 1, 3); // d(bus Ii eq)/d(fault Ii)

        return success.report(__func__);
      }

    private:
      /**
       * Checks the Jacobian value at a row and column is zero, treating an
       * absent entry as zero
       */
      bool jacobianEntryIsZero(const std::vector<DependencyTracking::Variable::DependencyMap>& jacobian,
                               size_t                                                          row,
                               size_t                                                          column) const
      {
        const auto entry = jacobian[row].find(column);
        if (entry == jacobian[row].end() || isEqual(entry->second, 0.0))
        {
          return true;
        }
        std::cout << "Jacobian entry (" << row << ", " << column
                  << ") = " << entry->second << " != 0\n";
        return false;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          const RealT R, const RealT X, const bool status)
      {
        DependencyTracking::Variable Vr1{1.0}; ///< Bus-1 real voltage
        DependencyTracking::Variable Vi1{1.0}; ///< Bus-1 imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>      bus(Vr1, Vi1);
        PhasorDynamics::BusFault<DependencyTracking::Variable, IdxT> fault(&bus, R, X, status);

        bus.allocate();
        fault.allocate();

        // Get d/dy
        bus.initialize();
        fault.initialize();

        auto* fault_y = fault.y().getData();
        for (size_t i = 0; i < fault.size(); ++i)
        {
          fault_y[i].setVariableNumber(i); ///< fault independent variables
        }
        fault.y().setDataUpdated();
        auto* bus_y = bus.y().getData();
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus_y[i].setVariableNumber(i + fault.size()); // Bus independent variables
        }
        bus.y().setDataUpdated();

        bus.evaluateResidual();
        fault.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                  ///< the dependencies
        auto&                                     residual_y_view = fault.getResidual();
        std::vector<DependencyTracking::Variable> residual_y(
            residual_y_view.getData(),
            residual_y_view.getData() + residual_y_view.getSize());

        // Get d/dy'
        bus.initialize();
        fault.initialize();

        auto* fault_yp = fault.yp().getData();
        for (size_t i = 0; i < fault.size(); ++i)
        {
          fault_yp[i].setVariableNumber(i); ///< fault independent variables
        }
        fault.yp().setDataUpdated();

        bus.evaluateResidual();
        fault.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                  ///< the dependencies
        auto&                                     residual_yp_view = fault.getResidual();
        std::vector<DependencyTracking::Variable> residual_yp(
            residual_yp_view.getData(),
            residual_yp_view.getData() + residual_yp_view.getSize());

        // Print the dependencies
        for (size_t i = 0; i < residual_y.size(); ++i)
        {
          std::cout << i << "th residual, y: ";
          (residual_y[i]).print(std::cout);
          std::cout << "\n";
          std::cout << i << "th residual, yp: ";
          (residual_yp[i]).print(std::cout);
          std::cout << "\n";
        }

        // Extract the dependencies and add d/dy' to d/dy
        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual_y.size());
        for (IdxT i = 0; i < residual_y.size(); ++i)
        {
          DependencyTracking::Variable::DependencyMap dependency_y  = (residual_y[i]).getDependencies();
          DependencyTracking::Variable::DependencyMap dependency_yp = (residual_yp[i]).getDependencies();

          for (const auto& pair_y : dependency_y)
          {
            auto index_y = pair_y.first;
            auto value_y = pair_y.second;
            auto it_yp   = dependency_yp.find(index_y);
            if (it_yp != dependency_yp.end())
            {
              auto value_yp = it_yp->second;
              dependencies[i].insert(std::make_pair(index_y, value_y + value_yp));
            }
            else
            {
              dependencies[i].insert(std::make_pair(index_y, value_y));
            }
          }

          // Insert yp dependencies that did not exist in the y dependencies
          for (const auto& pair_yp : dependency_yp)
          {
            auto index_yp = pair_yp.first;
            auto value_yp = pair_yp.second;
            auto it_y     = dependency_y.find(index_yp);
            if (it_y == dependency_y.end())
            {
              dependencies[i].insert(std::make_pair(index_yp, value_yp));
            }
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

        bus.initialize();
        fault.initialize();

        fault.updateTime(0.0, 1.0);

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + fault.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + fault.size()); // Reset bus residual indices
        }

        bus.evaluateResidual();
        fault.evaluateResidual();

        bus.evaluateJacobian();
        fault.evaluateJacobian();
        fault.constructCsr();
        GridKit::LinearAlgebra::CsrMatrix<ScalarT, IdxT>* model_jacobian = fault.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: BusFault Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif

    }; // class BusFaultTests

  } // namespace Testing
} // namespace GridKit
