#pragma once

#include <array>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class ExciterIeeet1Tests
    {
    public:
      using RealT                   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      static constexpr ScalarT tol_ = 100 * std::numeric_limits<ScalarT>::epsilon();

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

      TestOutcome zeroTimeConstantsAndDisabledSaturation()
      {
        TestStatus success = true;

        using Params = PhasorDynamics::Exciter::Ieeet1Parameters;

        auto data                    = makeTestData();
        data.parameters[Params::Tr]  = 0.0;
        data.parameters[Params::Ta]  = 0.0;
        data.parameters[Params::Ke]  = 0.1;
        data.parameters[Params::Te]  = 0.0;
        data.parameters[Params::Kf]  = 0.0;
        data.parameters[Params::Tf]  = 0.0;
        data.parameters[Params::E1]  = 0.0;
        data.parameters[Params::E2]  = 0.0;
        data.parameters[Params::Se1] = 0.0;
        data.parameters[Params::Se2] = 0.0;

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        ScalarT                                        efd_value{0.0};
        IdxT                                           efd_index = INVALID_INDEX<IdxT>;

        efd_node.set(&efd_value, &efd_index);
        exciter.getSignals()
            .template assignSignalNode<PhasorDynamics::Exciter::Ieeet1InternalVariables::EFD>(&efd_node);

        bus.allocate();
        exciter.allocate();

        bus.initialize();
        efd_node.init(1.2);
        success *= (exciter.initialize() == 0);
        exciter.tagDifferentiable();

        success *= (exciter.tag()[0]);
        success *= (exciter.tag()[1]);
        success *= (exciter.tag()[2]);
        success *= (exciter.tag()[3]);

        auto*       y  = exciter.y().getData();
        auto*       yp = exciter.yp().getData();
        const auto* f  = exciter.getResidual().getData();

        success *= isEqual(y[2], static_cast<ScalarT>(1.2), tol_);
        success *= isEqual(y[1], static_cast<ScalarT>(0.12), tol_);
        success *= isEqual(y[6], static_cast<ScalarT>(0.0), tol_);
        success *= isEqual(y[7], static_cast<ScalarT>(1.2), tol_);
        success *= isEqual(y[8], static_cast<ScalarT>(0.0), tol_);

        for (IdxT i = 0; i < exciter.y().getSize(); ++i)
        {
          success *= std::isfinite(y[i]);
        }

        exciter.evaluateResidual();
        for (IdxT i = 0; i < exciter.getResidual().getSize(); ++i)
        {
          success *= std::isfinite(f[i]);
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), tol_);
        }

        y[2] = 4.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[8], static_cast<ScalarT>(0.0), tol_);
        y[2]     = 1.2;

        yp[0] = 123.0;
        exciter.y().setDataUpdated();
        exciter.yp().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[0], static_cast<ScalarT>(-123.0), tol_);
        yp[0]    = 0.0;

        y[0] = 4.0;
        exciter.y().setDataUpdated();
        exciter.yp().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[0], static_cast<ScalarT>(1.0e3), tol_);

        y[0] = 5.0;
        y[4] = 0.02;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[1], static_cast<ScalarT>(880.0), tol_);

        y[4] = 0.0;
        y[1] = 1.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[2], static_cast<ScalarT>(880.0), tol_);

        y[1] = 0.0;
        y[5] = 1.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[3], static_cast<ScalarT>(1.0e3), tol_);

        return success.report(__func__);
      }

      TestOutcome automaticKeAndSpeedSelector()
      {
        TestStatus success = true;

        using Params = PhasorDynamics::Exciter::Ieeet1Parameters;

        auto run_case = [&](auto selector)
        {
          auto data                        = makeTestData();
          data.parameters[Params::Ke]      = 0.0;
          data.parameters[Params::Ispdlim] = selector;

          PhasorDynamics::Bus<ScalarT, IdxT>             bus(1.0, 0.0);
          PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
          PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT>      omega_node;
          ScalarT                                        efd_value{0.0};
          ScalarT                                        omega_value{0.1};
          IdxT                                           efd_index   = INVALID_INDEX<IdxT>;
          IdxT                                           omega_index = INVALID_INDEX<IdxT>;

          efd_node.set(&efd_value, &efd_index);
          omega_node.set(&omega_value, &omega_index);
          exciter.getSignals()
              .template assignSignalNode<PhasorDynamics::Exciter::Ieeet1InternalVariables::EFD>(&efd_node);
          exciter.getSignals()
              .template attachSignalNode<PhasorDynamics::Exciter::Ieeet1ExternalVariables::OMEGA>(&omega_node);

          bus.allocate();
          exciter.allocate();
          bus.initialize();
          efd_node.init(3.3);

          bool passed  = exciter.initialize() == 0;
          passed      &= exciter.evaluateResidual() == 0;

          const auto* y  = exciter.y().getData();
          passed        &= isEqual(y[2], static_cast<ScalarT>(3.0), tol_);
          passed        &= isEqual(y[1], static_cast<ScalarT>(0.1), tol_);
          passed        &= y[8] > static_cast<ScalarT>(0.0);

          const auto* f = exciter.getResidual().getData();
          for (IdxT i = 0; i < exciter.getResidual().getSize(); ++i)
          {
            passed &= isEqual(f[i], static_cast<ScalarT>(0.0), tol_);
          }
          return passed;
        };

        success *= run_case(true);
        success *= run_case(static_cast<IdxT>(2));
        success *= run_case(static_cast<RealT>(-1.0));

        auto zero_data                   = makeTestData();
        zero_data.parameters[Params::Ke] = 0.0;
        PhasorDynamics::Bus<ScalarT, IdxT>             zero_bus(1.0, 0.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> zero_exciter(&zero_bus, zero_data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      zero_efd_node;
        ScalarT                                        zero_efd{0.0};
        IdxT                                           zero_efd_index = INVALID_INDEX<IdxT>;
        zero_efd_node.set(&zero_efd, &zero_efd_index);
        zero_exciter.getSignals()
            .template assignSignalNode<PhasorDynamics::Exciter::Ieeet1InternalVariables::EFD>(&zero_efd_node);
        zero_bus.allocate();
        zero_exciter.allocate();
        zero_bus.initialize();
        success *= zero_exciter.initialize() != 0;

        return success.report(__func__);
      }

      TestOutcome saturationParameters()
      {
        TestStatus success = true;

        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing one-zero saturation fits and invalid saturation parameters. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);

        using Params = PhasorDynamics::Exciter::Ieeet1Parameters;

        struct FitCase
        {
          const char* label;
          RealT       e1;
          RealT       se1;
          RealT       e2;
          RealT       se2;
          RealT       knee;
          RealT       point;
          RealT       expected;
        };

        const std::array<FitCase, 2> fit_cases{{
            {"ascending one-zero fit", 2.8, 0.0, 3.373, 0.33, 2.8, 3.373, 0.33},
            {"descending one-zero fit", 3.373, 0.33, 2.8, 0.0, 2.8, 3.373, 0.33},
        }};

        for (const auto& test_case : fit_cases)
        {
          auto data                    = makeTestData();
          data.parameters[Params::E1]  = test_case.e1;
          data.parameters[Params::Se1] = test_case.se1;
          data.parameters[Params::E2]  = test_case.e2;
          data.parameters[Params::Se2] = test_case.se2;

          PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
          PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);

          bus.allocate();
          exciter.allocate();
          success *= (exciter.verify() == 0);

          auto* y = exciter.y().getData();
          y[2]    = test_case.knee;
          y[8]    = 0.0;
          exciter.y().setDataUpdated();
          exciter.evaluateResidual();
          success *= isEqual(exciter.getResidual().getData()[8], static_cast<ScalarT>(0.0));

          y[2] = test_case.point;
          y[8] = 0.0;
          exciter.y().setDataUpdated();
          exciter.evaluateResidual();
          if (!isEqual(exciter.getResidual().getData()[8],
                       static_cast<ScalarT>(test_case.expected)))
          {
            std::cout << test_case.label << " did not reproduce its nonzero point\n";
            success = false;
          }
        }

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);

        auto negative                    = makeTestData();
        negative.parameters[Params::Se1] = -0.1;
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> negative_exciter(&bus, negative);
        bus.allocate();
        negative_exciter.allocate();
        auto* negative_y   = negative_exciter.y().getData();
        auto* negative_yp  = negative_exciter.yp().getData();
        negative_y[0]      = static_cast<ScalarT>(17.0);
        negative_yp[0]     = static_cast<ScalarT>(19.0);
        success           *= (negative_exciter.verify() != 0);
        success           *= (negative_exciter.initialize() != 0);
        success           *= (negative_y[0] == static_cast<ScalarT>(17.0));
        success           *= (negative_yp[0] == static_cast<ScalarT>(19.0));

        auto crossed                    = makeTestData();
        crossed.parameters[Params::E1]  = 3.373;
        crossed.parameters[Params::Se1] = 0.0;
        crossed.parameters[Params::E2]  = 2.8;
        crossed.parameters[Params::Se2] = 0.33;
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> crossed_exciter(&bus, crossed);
        success *= (crossed_exciter.verify() != 0);

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
        exciter.y().setDataUpdated();
        auto* bus_y = bus.y().getData();
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus_y[i].setVariableNumber(i + exciter.size()); // Bus independent variables
        }
        bus.y().setDataUpdated();

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
        exciter.yp().setDataUpdated();

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
