#pragma once

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class GovernorGastPtiTests
    {
    public:
      using RealT  = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using Gov    = PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;
      using Data   = PhasorDynamics::Governor::GastPtiData<RealT, IdxT>;
      using Var    = PhasorDynamics::Governor::GastPtiInternalVariables;
      using Ext    = PhasorDynamics::Governor::GastPtiExternalVariables;
      using Params = PhasorDynamics::Governor::GastPtiParameters;
      using Mon    = PhasorDynamics::Governor::GastPtiMonitorableVariables;

      GovernorGastPtiTests()  = default;
      ~GovernorGastPtiTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        Gov model(makeTestData());

        const auto* monitor  = model.getMonitor();
        success             *= (model.size() == static_cast<IdxT>(Var::MAXIMUM));
        success             *= (monitor != nullptr);
        if (monitor != nullptr)
        {
          success *= (!monitor->empty());
        }
        success *= (model.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        Gov model(makeTestData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        const ScalarT                             pmech0 = scalar(kInitialPmech);
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   omega_value = scalar(kInitialOmega);
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      omega_index = 9;
        pmech_node.set(&pmech_value, &pmech_index);
        omega_node.set(&omega_value, &omega_index);

        model.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);
        model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);

        success *= (model.allocate() == 0);
        pmech_node.init(pmech0);
        success *= pmech_node.linked();
        success *= (pmech_node.getVariableIndex() == static_cast<IdxT>(Var::PMECH));

        success *= (model.verify() == 0);
        success *= (model.initialize() == 0);
        success *= (model.tagDifferentiable() == 0);
        success *= (model.evaluateResidual() == 0);

        const ScalarT xflow0 = pmech0 + scalar(kDturb) * omega_value;
        const ScalarT vtemp0 = scalar(kAt) + scalar(kKt) * (scalar(kAt) - xflow0);

        success *= isEqual(model.y()[index(Var::XVALVE)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y()[index(Var::XFLOW)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y()[index(Var::XTEMP)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y()[index(Var::VLOAD)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y()[index(Var::VTEMP)], vtemp0, scalar(kTolerance));
        success *= isEqual(model.y()[index(Var::VLV)], xflow0, scalar(kTolerance));
        success *= (model.tag()[index(Var::XVALVE)] == true);
        success *= (model.tag()[index(Var::XFLOW)] == true);
        success *= (model.tag()[index(Var::XTEMP)] == true);

        checkZeroResidual(model, success);

        return success.report(__func__);
      }

      TestOutcome prefSignal()
      {
        TestStatus success = true;

        Gov model(makeTestData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        const ScalarT                             pmech0 = scalar(kInitialPmech);
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   omega_value = scalar(kInitialOmega);
        ScalarT                                   pref_value  = prefForInitialPoint(pmech0, omega_value);
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      omega_index = 9;
        IdxT                                      pref_index  = 10;
        pmech_node.set(&pmech_value, &pmech_index);
        omega_node.set(&omega_value, &omega_index);
        pref_node.set(&pref_value, &pref_index);

        model.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);
        model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);

        success *= (model.allocate() == 0);
        pmech_node.init(pmech0);
        success *= (model.verify() == 0);
        success *= (model.initialize() == 0);
        success *= (model.evaluateResidual() == 0);
        checkZeroResidual(model, success);

        pref_value += scalar(kPrefStep);
        success    *= (model.evaluateResidual() == 0);
        success    *= isEqual(model.getResidual()[index(Var::VLOAD)],
                           scalar(kPrefStep),
                           scalar(kTolerance));

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        Gov model(makeTestData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        ScalarT                                   omega_value = scalar(kResidualOmega);
        ScalarT                                   pref_value  = scalar(kResidualPref);
        IdxT                                      omega_index = 7;
        IdxT                                      pref_index  = 8;
        omega_node.set(&omega_value, &omega_index);
        pref_node.set(&pref_value, &pref_index);

        model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);

        success *= (model.allocate() == 0);

        const ScalarT xvalve     = scalar(kResidualXvalve);
        const ScalarT xflow      = scalar(kResidualXflow);
        const ScalarT xtemp      = scalar(kResidualXtemp);
        const ScalarT vload      = scalar(kResidualVload);
        const ScalarT vtemp      = scalar(kResidualVtemp);
        const ScalarT vlv        = scalar(kResidualVlv);
        const ScalarT pmech      = scalar(kResidualPmech);
        const ScalarT xvalve_dot = scalar(kResidualXvalveDot);
        const ScalarT xflow_dot  = scalar(kResidualXflowDot);
        const ScalarT xtemp_dot  = scalar(kResidualXtempDot);

        model.y()[index(Var::XVALVE)] = xvalve;
        model.y()[index(Var::XFLOW)]  = xflow;
        model.y()[index(Var::XTEMP)]  = xtemp;
        model.y()[index(Var::VLOAD)]  = vload;
        model.y()[index(Var::VTEMP)]  = vtemp;
        model.y()[index(Var::VLV)]    = vlv;
        model.y()[index(Var::PMECH)]  = pmech;

        model.yp()[index(Var::XVALVE)] = xvalve_dot;
        model.yp()[index(Var::XFLOW)]  = xflow_dot;
        model.yp()[index(Var::XTEMP)]  = xtemp_dot;

        success *= (model.verify() == 0);
        success *= (model.evaluateResidual() == 0);

        const ScalarT              valve_target = vlv - xvalve;
        const ScalarT              selected_vlv = vload;
        const std::vector<ScalarT> expected     = {
            -scalar(kT1) * xvalve_dot + valve_target,
            -scalar(kT2) * xflow_dot - xflow + xvalve,
            -scalar(kT3) * xtemp_dot - xtemp + xflow,
            -vload + pref_value - omega_value / scalar(kR),
            -vtemp + scalar(kAt) + scalar(kKt) * (scalar(kAt) - xtemp),
            -vlv + selected_vlv,
            -pmech + xflow - scalar(kDturb) * omega_value,
        };

        checkResidual(model, expected, success);

        return success.report(__func__);
      }

      TestOutcome antiWindupLimiter()
      {
        TestStatus success = true;

        Gov model(makeTestData());
        success *= (model.allocate() == 0);

        auto check_valve = [&](ScalarT xvalve, ScalarT vlv, ScalarT expected)
        {
          std::fill(model.y().begin(), model.y().end(), ScalarT{0});
          std::fill(model.yp().begin(), model.yp().end(), ScalarT{0});
          model.y()[index(Var::XVALVE)]  = xvalve;
          model.y()[index(Var::VLV)]     = vlv;
          success                       *= (model.evaluateResidual() == 0);
          success                       *= isEqual(model.getResidual()[index(Var::XVALVE)],
                             expected,
                             scalar(kSmoothTolerance));
        };

        check_valve(scalar(2.2), scalar(3.2), scalar(0.0));
        check_valve(scalar(2.2), scalar(1.2), scalar(-1.0));
        check_valve(scalar(-1.0), scalar(-2.0), scalar(0.0));
        check_valve(scalar(-1.0), scalar(0.0), scalar(1.0));

        return success.report(__func__);
      }

      TestOutcome timeConstantTags()
      {
        TestStatus success = true;

        auto data                   = makeTestData();
        data.parameters[Params::T1] = static_cast<RealT>(0.0);
        data.parameters[Params::T2] = static_cast<RealT>(kT2);
        data.parameters[Params::T3] = static_cast<RealT>(0.0);

        Gov model(data);
        success *= (model.allocate() == 0);
        success *= (model.tagDifferentiable() == 0);
        success *= (model.tag()[index(Var::XVALVE)] == false);
        success *= (model.tag()[index(Var::XFLOW)] == true);
        success *= (model.tag()[index(Var::XTEMP)] == false);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        auto missing = makeTestData();
        missing.parameters.erase(Params::R);
        Gov missing_model(missing);
        success *= (missing_model.verify() > 0);

        auto negative_time                   = makeTestData();
        negative_time.parameters[Params::T2] = static_cast<RealT>(-kT2);
        Gov negative_time_model(negative_time);
        success *= (negative_time_model.verify() > 0);

        auto invalid_limits                     = makeTestData();
        invalid_limits.parameters[Params::Vmin] = static_cast<RealT>(kVmax + 0.1);
        invalid_limits.parameters[Params::Vmax] = static_cast<RealT>(kVmax);
        Gov invalid_limits_model(invalid_limits);
        success *= (invalid_limits_model.verify() > 0);

        auto invalid_optional                      = makeTestData();
        invalid_optional.parameters[Params::Trate] = true;
        Gov invalid_optional_model(invalid_optional);
        success *= (invalid_optional_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome signalValidation()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        Gov                                       omega_model(makeTestData());
        omega_model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        success *= (omega_model.verify() > 0);

        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        Gov                                       pref_model(makeTestData());
        pref_model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        success *= (pref_model.verify() > 0);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto dependency_tracking_jacobian = dependencyTrackingJacobian();
        auto enzyme_jacobian              = enzymeJacobian();

        success *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], kTolerance);
        }

        return success.report(__func__);
      }
#endif

    private:
      static constexpr RealT kTolerance       = 100.0 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT kSmoothTolerance = 1.0e-2;

      static constexpr RealT kR     = 0.05;
      static constexpr RealT kT1    = 0.4;
      static constexpr RealT kT2    = 0.5;
      static constexpr RealT kT3    = 0.25;
      static constexpr RealT kAt    = 2.0;
      static constexpr RealT kKt    = 0.3;
      static constexpr RealT kVmax  = 1.2;
      static constexpr RealT kVmin  = 0.0;
      static constexpr RealT kDturb = 0.1;
      static constexpr RealT kTrate = 0.0;

      static constexpr RealT kInitialPmech = 0.75;
      static constexpr RealT kInitialOmega = 0.02;
      static constexpr RealT kPrefStep     = 0.1;

      static constexpr RealT kResidualOmega     = 0.02;
      static constexpr RealT kResidualPref      = 1.25;
      static constexpr RealT kResidualXvalve    = 0.7;
      static constexpr RealT kResidualXflow     = 0.6;
      static constexpr RealT kResidualXtemp     = 0.5;
      static constexpr RealT kResidualVload     = 0.9;
      static constexpr RealT kResidualVtemp     = 2.5;
      static constexpr RealT kResidualVlv       = 0.8;
      static constexpr RealT kResidualPmech     = 0.55;
      static constexpr RealT kResidualXvalveDot = 0.05;
      static constexpr RealT kResidualXflowDot  = -0.1;
      static constexpr RealT kResidualXtempDot  = 0.2;

      static ScalarT scalar(RealT value)
      {
        return static_cast<ScalarT>(value);
      }

      template <class ValueT>
      static ValueT value(RealT value)
      {
        return ValueT{value};
      }

      static size_t index(Var variable)
      {
        return static_cast<size_t>(variable);
      }

      template <class ValueT>
      static ValueT prefForInitialPoint(const ValueT& pmech, const ValueT& omega)
      {
        return pmech + value<ValueT>(kDturb) * omega + omega / value<ValueT>(kR);
      }

      void checkResidual(const Gov&                  model,
                         const std::vector<ScalarT>& expected,
                         TestStatus&                 success) const
      {
        const auto& residual = model.getResidual();
        for (size_t i = 0; i < expected.size(); ++i)
        {
          if (!isEqual(residual[i], expected[i], scalar(kTolerance)))
          {
            std::cout << "Unexpected GASTPTI residual at index " << i << ": "
                      << residual[i] << " != " << expected[i] << "\n";
            success = false;
          }
        }
      }

      void checkZeroResidual(const Gov& model, TestStatus& success) const
      {
        std::vector<ScalarT> expected(static_cast<size_t>(Var::MAXIMUM), ScalarT{0});
        checkResidual(model, expected, success);
      }

      Data makeTestData()
      {
        Data data;
        data.device_class          = "GastPti";
        data.disambiguation_string = "gastpti_test";
        data.monitored_variables.insert(Mon::pmech);
        data.monitored_variables.insert(Mon::fuelvalve);

        data.parameters[Params::R]     = static_cast<RealT>(kR);
        data.parameters[Params::T1]    = static_cast<RealT>(kT1);
        data.parameters[Params::T2]    = static_cast<RealT>(kT2);
        data.parameters[Params::T3]    = static_cast<RealT>(kT3);
        data.parameters[Params::At]    = static_cast<RealT>(kAt);
        data.parameters[Params::Kt]    = static_cast<RealT>(kKt);
        data.parameters[Params::Vmax]  = static_cast<RealT>(kVmax);
        data.parameters[Params::Vmin]  = static_cast<RealT>(kVmin);
        data.parameters[Params::Dturb] = static_cast<RealT>(kDturb);
        data.parameters[Params::Trate] = static_cast<RealT>(kTrate);

        return data;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      using DependencyMap = DependencyTracking::Variable::DependencyMap;

      std::vector<DependencyMap> dependencyTrackingJacobian()
      {
        using ADScalarT = DependencyTracking::Variable;
        using ADGov     = PhasorDynamics::Governor::GastPti<ADScalarT, IdxT>;

        ADGov model(makeTestData());

        PhasorDynamics::SignalNode<ADScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ADScalarT, IdxT> pref_node;
        ADScalarT                                   omega_value{kResidualOmega};
        ADScalarT                                   pref_value{kResidualPref};
        IdxT                                        omega_index = static_cast<IdxT>(Var::MAXIMUM);
        IdxT                                        pref_index  = omega_index + 1;
        omega_node.set(&omega_value, &omega_index);
        pref_node.set(&pref_value, &pref_index);

        model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        model.allocate();
        model.updateTime(0.0, 1.0);

        model.y()[index(Var::XVALVE)] = ADScalarT{kResidualXvalve};
        model.y()[index(Var::XFLOW)]  = ADScalarT{kResidualXflow};
        model.y()[index(Var::XTEMP)]  = ADScalarT{kResidualXtemp};
        model.y()[index(Var::VLOAD)]  = ADScalarT{kResidualVload};
        model.y()[index(Var::VTEMP)]  = ADScalarT{kResidualVtemp};
        model.y()[index(Var::VLV)]    = ADScalarT{kResidualVlv};
        model.y()[index(Var::PMECH)]  = ADScalarT{kResidualPmech};

        model.yp()[index(Var::XVALVE)] = ADScalarT{kResidualXvalveDot};
        model.yp()[index(Var::XFLOW)]  = ADScalarT{kResidualXflowDot};
        model.yp()[index(Var::XTEMP)]  = ADScalarT{kResidualXtempDot};

        for (size_t i = 0; i < model.y().size(); ++i)
        {
          model.y()[i].setVariableNumber(i);
          model.yp()[i].setVariableNumber(i);
        }
        omega_value.setVariableNumber(static_cast<size_t>(omega_index));
        pref_value.setVariableNumber(static_cast<size_t>(pref_index));

        model.evaluateResidual();

        std::vector<DependencyMap> dependencies(model.getResidual().size());
        for (size_t i = 0; i < dependencies.size(); ++i)
        {
          dependencies[i] = model.getResidual()[i].getDependencies();
        }

        return dependencies;
      }

      std::vector<DependencyMap> enzymeJacobian()
      {
        Gov model(makeTestData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        ScalarT                                   omega_value = scalar(kResidualOmega);
        ScalarT                                   pref_value  = scalar(kResidualPref);
        IdxT                                      omega_index = static_cast<IdxT>(Var::MAXIMUM);
        IdxT                                      pref_index  = omega_index + 1;
        omega_node.set(&omega_value, &omega_index);
        pref_node.set(&pref_value, &pref_index);

        model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        model.allocate();
        model.updateTime(0.0, 1.0);

        model.y()[index(Var::XVALVE)] = scalar(kResidualXvalve);
        model.y()[index(Var::XFLOW)]  = scalar(kResidualXflow);
        model.y()[index(Var::XTEMP)]  = scalar(kResidualXtemp);
        model.y()[index(Var::VLOAD)]  = scalar(kResidualVload);
        model.y()[index(Var::VTEMP)]  = scalar(kResidualVtemp);
        model.y()[index(Var::VLV)]    = scalar(kResidualVlv);
        model.y()[index(Var::PMECH)]  = scalar(kResidualPmech);

        model.yp()[index(Var::XVALVE)] = scalar(kResidualXvalveDot);
        model.yp()[index(Var::XFLOW)]  = scalar(kResidualXflowDot);
        model.yp()[index(Var::XTEMP)]  = scalar(kResidualXtempDot);

        model.evaluateResidual();
        model.evaluateJacobian();

        auto& model_jacobian = model.getJacobian();
        model_jacobian.deduplicate();

        return MapFromCOO(model_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
