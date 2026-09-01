#pragma once

#include <iostream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
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
      using RealT    = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using External = PhasorDynamics::BusFaultExternalVariables;

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
       * Verifies the fault current injected into the bus residual
       */
      TestOutcome residual(bool status)
      {
        TestStatus success = true;

        const RealT R{0.1};
        const RealT X{1e-3};

        const ScalarT Vr{1.0}; ///< Bus real voltage
        const ScalarT Vi{0.5}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr, Vi);
        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        PhasorDynamics::BusFaultData<RealT, IdxT> fault_data;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::R] = R;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::X] = X;
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, fault_data);

        ScalarT status_value{0.0};
        if (status)
        {
          status_value = 1.0;
        }
        IdxT                                      status_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<ScalarT, IdxT> status_node;
        status_node.set(&status_value, &status_index);
        fault.getSignals().template attachSignalNode<External::STATUS>(&status_node);

        fault.allocate();
        fault.evaluateResidual();

        const auto [g, b] = faultAdmittance(R, X, status);

        const ScalarT Ir_expected{g * Vr - b * Vi};
        const ScalarT Ii_expected{b * Vr + g * Vi};

        success *= isEqual(bus.Ir(), Ir_expected);
        success *= isEqual(bus.Ii(), Ii_expected);

        return success.report(__func__);
      }

      /**
       * A test case to verify Jacobian values via dependency tracking
       */
      TestOutcome jacobian(bool status)
      {
        TestStatus success = true;

        const RealT R{0.1};
        const RealT X{1e-3};

        DependencyTracking::Variable Vr{1.0}; ///< Bus real voltage
        DependencyTracking::Variable Vi{0.5}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT> bus(Vr, Vi);
        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();
        bus.Vr().setVariableNumber(0);
        bus.Vi().setVariableNumber(1);

        PhasorDynamics::BusFaultData<RealT, IdxT> fault_data;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::R] = R;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::X] = X;
        PhasorDynamics::BusFault<DependencyTracking::Variable, IdxT> fault(&bus, fault_data);

        DependencyTracking::Variable status_value{0.0};
        if (status)
        {
          status_value = DependencyTracking::Variable{1.0};
        }
        IdxT                                                           status_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT> status_node;
        status_node.set(&status_value, &status_index);
        fault.getSignals().template attachSignalNode<External::STATUS>(&status_node);

        fault.allocate();
        fault.evaluateResidual(); ///< Computes the residual and the Jacobian values
                                  ///< by tracking the dependencies

        std::vector<DependencyTracking::Variable>                residuals{bus.Ir(), bus.Ii()};
        std::vector<DependencyTracking::Variable::DependencyMap> ref =
            analyticalJacobian(R, X, status);

        /// Compare dependencies computed automatically to the ones computed
        /// analytically
        for (size_t i = 0; i < residuals.size(); ++i)
        {
          const DependencyTracking::Variable::DependencyMap& dependencies =
              residuals[i].getDependencies();
          success *= (GridKit::Testing::isEqual(dependencies, ref[i]));
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * A test case to verify Enzyme Jacobian values against the analytical
       * Jacobian
       */
      TestOutcome enzymeJacobian(bool status)
      {
        TestStatus success = true;

        const RealT R{0.1};
        const RealT X{1e-3};

        const ScalarT Vr{1.0}; ///< Bus real voltage
        const ScalarT Vi{0.5}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr, Vi);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::BusFaultData<RealT, IdxT> fault_data;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::R] = R;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::X] = X;
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, fault_data);

        ScalarT status_value{0.0};
        if (status)
        {
          status_value = 1.0;
        }
        IdxT                                      status_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<ScalarT, IdxT> status_node;
        status_node.set(&status_value, &status_index);
        fault.getSignals().template attachSignalNode<External::STATUS>(&status_node);

        fault.allocate();

        bus.evaluateResidual();
        fault.evaluateResidual();

        bus.evaluateJacobian();
        fault.evaluateJacobian();
        fault.constructCsr();
        GridKit::LinearAlgebra::CsrMatrix<ScalarT, IdxT>* model_jacobian =
            fault.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: BusFault Jacobian\n";
        model_jacobian->print();

        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian =
            GridKit::Testing::MapFromCsr(model_jacobian);
        std::vector<DependencyTracking::Variable::DependencyMap> ref =
            analyticalJacobian(R, X, status);

        /// Compare Enzyme dependencies to the ones computed analytically
        for (size_t i = 0; i < ref.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(enzyme_jacobian[i], ref[i]));
        }

        return success.report(__func__);
      }

      /**
       * Verify the fixed Jacobian pattern and status-scaled values across a
       * complete fault event.
       */
      TestOutcome jacobianAcrossStatusChanges()
      {
        TestStatus success = true;

        const RealT R{0.0};
        const RealT X{0.01};

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::BusFaultData<RealT, IdxT> fault_data;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::R] = R;
        fault_data.parameters[PhasorDynamics::BusFaultParameters::X] = X;
        PhasorDynamics::BusFault<ScalarT, IdxT> fault(&bus, fault_data);

        ScalarT                                   status_value{0.0};
        IdxT                                      status_index{INVALID_INDEX<IdxT>};
        PhasorDynamics::SignalNode<ScalarT, IdxT> status_node;
        status_node.set(&status_value, &status_index);
        fault.getSignals().template attachSignalNode<External::STATUS>(&status_node);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        system.addBus(&bus);
        system.addFault(&fault);
        system.allocate();
        system.initialize();

        auto checkJacobian = [&](const bool status)
        {
          system.evaluateResidual();
          system.evaluateJacobian();

          const auto jacobian = GridKit::Testing::MapFromCsr(system.getCsrJacobian());
          const auto [g, b]   = faultAdmittance(R, X, status);

          success *= jacobianEntryEquals(jacobian, 0, 0, g);
          success *= jacobianEntryEquals(jacobian, 0, 1, -b);
          success *= jacobianEntryEquals(jacobian, 1, 0, b);
          success *= jacobianEntryEquals(jacobian, 1, 1, g);
        };

        checkJacobian(false);

        status_node.init(1.0);
        checkJacobian(true);

        status_node.init(0.0);
        checkJacobian(false);

        return success.report(__func__);
      }
#endif

    private:
      /**
       * Checks that a Jacobian entry exists and has the expected value.
       */
      bool jacobianEntryEquals(const std::vector<DependencyTracking::Variable::DependencyMap>& jacobian,
                               const size_t                                                    row,
                               const size_t                                                    column,
                               const RealT                                                     expected) const
      {
        if (row >= jacobian.size())
        {
          std::cout << "Jacobian row " << row << " is missing\n";
          return false;
        }

        const auto entry = jacobian[row].find(column);
        if (entry == jacobian[row].end())
        {
          std::cout << "Jacobian entry (" << row << ", " << column << ") is missing\n";
          return false;
        }

        if (!isEqual(entry->second, expected))
        {
          std::cout << "Jacobian entry (" << row << ", " << column
                    << ") = " << entry->second << " != " << expected << "\n";
          return false;
        }

        return true;
      }

      std::pair<RealT, RealT> faultAdmittance(const RealT R, const RealT X, const bool status)
      {
        RealT g{0.0};
        RealT b{0.0};

        if (status)
        {
          const RealT denom = R * R + X * X;
          g                 = -R / denom;
          b                 = X / denom;
        }

        return {g, b};
      }

      std::vector<DependencyTracking::Variable::DependencyMap>
      analyticalJacobian(const RealT R, const RealT X, const bool status)
      {
        const auto [g, b] = faultAdmittance(R, X, status);

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(2);
        dependencies[0] = {{0, g}, {1, -b}};
        dependencies[1] = {{0, b}, {1, g}};

        return dependencies;
      }

    }; // class BusFaultTests

  } // namespace Testing
} // namespace GridKit
