#pragma once

#include <cmath>
#include <iostream>
#include <limits>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ExciterIeeet1Tests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ExciterIeeet1Tests()  = default;
      ~ExciterIeeet1Tests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        auto                               data    = makeTestData();
        auto*                              exciter = new PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT>(&bus, data);

        success *= (exciter != nullptr);
        success *= (exciter->size() == 9);
        success *= (exciter->getMonitor() != nullptr);

        delete exciter;

        return success.report(__func__);
      }

      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);

        bus.allocate();
        exciter.allocate();

        bus.initialize();
        exciter.initialize();

        exciter.evaluateResidual();

        const auto& residual      = exciter.getResidual();
        const auto* residual_data = residual.getData();
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          if (!isEqual(residual_data[i], static_cast<ScalarT>(0.0)))
          {
            std::cout << "Non-zero Ieeet1 residual at index " << i << ": " << residual_data[i] << "\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

      TestOutcome zeroTrUsesMinimum()
      {
        TestStatus success = true;

        auto data                                                      = makeTestData();
        data.parameters[PhasorDynamics::Exciter::Ieeet1Parameters::Tr] = 0.0;

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);

        bus.allocate();
        exciter.allocate();

        bus.initialize();
        exciter.initialize();
        exciter.tagDifferentiable();

        const auto* tag  = exciter.tag().getData();
        auto*       y    = exciter.y().getData();
        auto*       yp   = exciter.yp().getData();
        success         *= tag[0];

        yp[0] = 123.0;
        exciter.evaluateResidual();
        const auto* f  = exciter.getResidual().getData();
        success       *= isEqual(f[0], static_cast<ScalarT>(-123.0));

        yp[0] = 0.0;
        y[0]  = 4.0;
        exciter.evaluateResidual();
        success *= isEqual(f[0], static_cast<ScalarT>(1.0e3));

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * @brief Checks Jacobian evaluation.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto tol = 10 * std::numeric_limits<RealT>::epsilon();

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian();

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian();

        // Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        auto data = makeTestData();

        DependencyTracking::Variable                                        Vr1{3.0};
        DependencyTracking::Variable                                        Vi1{4.0};
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>             bus(Vr1, Vi1);
        PhasorDynamics::Exciter::Ieeet1<DependencyTracking::Variable, IdxT> exciter(&bus, data);

        bus.allocate();
        exciter.allocate();

        // Get d/dy
        bus.initialize();
        exciter.initialize();

        auto* exciter_y = exciter.y().getData();
        for (size_t i = 0; i < exciter.size(); ++i)
        {
          exciter_y[i].setVariableNumber(i); ///< Exciter independent variables
        }
        auto* bus_y = bus.y().getData();
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus_y[i].setVariableNumber(i + exciter.size()); // Bus independent variables
        }

        bus.evaluateResidual();
        exciter.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                    ///< the dependencies
        auto&                                     residual_y_view = exciter.getResidual();
        std::vector<DependencyTracking::Variable> residual_y(residual_y_view.getData(), residual_y_view.getData() + residual_y_view.getSize());

        // Get d/dy'
        bus.initialize();
        exciter.initialize();

        auto* exciter_yp = exciter.yp().getData();
        for (size_t i = 0; i < exciter.size(); ++i)
        {
          exciter_yp[i].setVariableNumber(i); ///< Exciter independent variables
        }

        bus.evaluateResidual();
        exciter.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                    ///< the dependencies
        auto&                                     residual_yp_view = exciter.getResidual();
        std::vector<DependencyTracking::Variable> residual_yp(residual_yp_view.getData(), residual_yp_view.getData() + residual_yp_view.getSize());

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

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);

        bus.allocate();
        exciter.allocate();

        bus.initialize();
        exciter.initialize();

        exciter.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + exciter.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + exciter.size()); // Reset bus residual indices
        }

        bus.evaluateResidual();
        exciter.evaluateResidual();

        bus.evaluateJacobian();
        exciter.evaluateJacobian();
        exciter.constructCsr();
        GridKit::LinearAlgebra::CsrMatrix<ScalarT, IdxT>* model_jacobian = exciter.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: Ieeet1 Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif

    private:
      auto makeTestData() -> PhasorDynamics::Exciter::Ieeet1Data<RealT, IdxT>
      {
        using Params = PhasorDynamics::Exciter::Ieeet1Parameters;

        PhasorDynamics::Exciter::Ieeet1Data<RealT, IdxT> data;
        data.device_class          = "exciter";
        data.disambiguation_string = "ieeet1_test";
        data.monitored_variables.insert(PhasorDynamics::Exciter::Ieeet1MonitorableVariables::efd);

        data.parameters[Params::Tr]      = 0.0;
        data.parameters[Params::Ka]      = 50.0;
        data.parameters[Params::Ta]      = 0.04;
        data.parameters[Params::Ke]      = -0.06;
        data.parameters[Params::Te]      = 0.6;
        data.parameters[Params::Kf]      = 0.09;
        data.parameters[Params::Tf]      = 1.46;
        data.parameters[Params::Vrmin]   = -1.0;
        data.parameters[Params::Vrmax]   = 1.0;
        data.parameters[Params::E1]      = 2.8;
        data.parameters[Params::E2]      = 3.373;
        data.parameters[Params::Se1]     = 0.04;
        data.parameters[Params::Se2]     = 0.33;
        data.parameters[Params::Ispdlim] = 0.0;

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
