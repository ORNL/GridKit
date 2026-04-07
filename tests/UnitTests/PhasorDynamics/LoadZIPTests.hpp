#pragma once

#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIP.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class LoadZIPTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      LoadZIPTests()  = default;
      ~LoadZIPTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* load =
            new PhasorDynamics::LoadZIP<ScalarT, IdxT>(bus);

        success *= (load != nullptr);

        if (load)
        {
          delete load;
        }
        delete bus;

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        RealT P0{2.0};
        RealT Q0{0.5};
        RealT V0{2.0};
        RealT alphaI{4.0};
        RealT alphaP{2.0};

        ScalarT Vr{0.3}; ///< Bus real voltage
        ScalarT Vi{0.4}; ///< Bus imaginary voltage

        const ScalarT Ir{-10.88}; ///< Solution real current
        const ScalarT Ii{-8.84};  ///< Solution imaginary current

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(Vr, Vi);
        PhasorDynamics::LoadZIP<ScalarT, IdxT>     load(&bus, P0, Q0, V0, alphaI, alphaP);
        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        load.evaluateResidual();

        success *= isEqual(bus.Ir(), Ir);
        success *= isEqual(bus.Ii(), Ii);

        success.expectFailure(); // TODO: Check the reference solution values
                                 // and the implementation

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * A test case to verify Jacobian values
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        RealT P0{2.0};
        RealT Q0{0.5};
        RealT V0{2.0};
        RealT alphaI{4.0};
        RealT alphaP{2.0};

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap>
            dependency_tracking_jacobian = DependencyTrackingJacobian(P0, Q0, V0, alphaI, alphaP);

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap>
            enzyme_jacobian = EnzymeJacobian(P0, Q0, V0, alphaI, alphaP);

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]));
        }

        success.expectFailure(); // TODO: Activate df/dy term in Enzyme Jacobian

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          const RealT P0, const RealT Q0, const RealT V0, const RealT alphaI, const RealT alphaP)
      {
        DependencyTracking::Variable Vr{0.3}; ///< Bus real voltage
        DependencyTracking::Variable Vi{0.4}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>     bus(Vr, Vi);
        PhasorDynamics::LoadZIP<DependencyTracking::Variable, IdxT> load(&bus, P0, Q0, V0, alphaI, alphaP);

        bus.allocate();
        load.allocate();

        // Get d/dy
        bus.initialize();
        load.initialize();

        for (size_t i = 0; i < load.size(); ++i)
        {
          load.y()[i].setVariableNumber(i); ///< load independent variables
        }
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.y()[i].setVariableNumber(i + load.size()); // Bus independent variables
        }

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies
        std::vector<DependencyTracking::Variable> residual_y = load.getResidual();

        // Get d/dy'
        bus.initialize();
        load.initialize();

        for (size_t i = 0; i < load.size(); ++i)
        {
          load.yp()[i].setVariableNumber(i); ///< load independent variables
        }

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies
        std::vector<DependencyTracking::Variable> residual_yp = load.getResidual();

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
          const RealT P0, const RealT Q0, const RealT V0, const RealT alphaI, const RealT alphaP)
      {
        ScalarT Vr{0.3}; ///< Bus real voltage
        ScalarT Vi{0.4}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>     bus(Vr, Vi);
        PhasorDynamics::LoadZIP<ScalarT, IdxT> load(&bus, P0, Q0, V0, alphaI, alphaP);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        load.updateTime(0.0, 1.0);

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + load.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + load.size()); // Reset bus residual indices
        }

        bus.evaluateResidual();
        load.evaluateResidual();

        bus.evaluateJacobian();
        load.evaluateJacobian();
        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& model_jacobian = load.getJacobian();
        model_jacobian.deduplicate();
        model_jacobian.printMatrix("Model Jacobian");

        return GridKit::Testing::MapFromCOO(model_jacobian);
      }
#endif
    }; // class LoadZIPTests

  } // namespace Testing
} // namespace GridKit
