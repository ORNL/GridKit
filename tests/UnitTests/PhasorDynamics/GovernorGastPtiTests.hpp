#pragma once

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class GovernorGastPtiTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using Gov     = PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;
      using Data    = PhasorDynamics::Governor::GastPtiData<RealT, IdxT>;
      using Var     = PhasorDynamics::Governor::GastPtiInternalVariables;
      using Ext     = PhasorDynamics::Governor::GastPtiExternalVariables;
      using Params  = PhasorDynamics::Governor::GastPtiParameters;
      using Mon     = PhasorDynamics::Governor::GastPtiMonitorableVariables;
      using Mode    = PhasorDynamics::Governor::ResponseMode;

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
        IdxT                                      omega_index = static_cast<IdxT>(Var::MAXIMUM);
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

        success *= isEqual(model.y().getData()[index(Var::XVALVE)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::XFLOW)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::XTEMP)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VLOAD)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VTEMP)], vtemp0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VLV)], xflow0, scalar(kTolerance));
        success *= (model.tag()[index(Var::XVALVE)] == true);
        success *= (model.tag()[index(Var::XFLOW)] == true);
        success *= (model.tag()[index(Var::XTEMP)] == true);

        checkZeroResidual(model, success);

        return success.report(__func__);
      }

      TestOutcome baseConversion()
      {
        TestStatus success = true;

        auto data                      = makeTestData();
        data.parameters[Params::Trate] = static_cast<RealT>(kConversionTrate);
        Gov model(data);
        model.setSystemBase(static_cast<RealT>(kSystemFrequency),
                            static_cast<RealT>(kConversionSystemBase * 1.0e6));

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        const ScalarT                             pmech0 = scalar(kConversionInitialPmech);
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   omega_value = scalar(kInitialOmega);
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      omega_index = static_cast<IdxT>(Var::MAXIMUM);
        pmech_node.set(&pmech_value, &pmech_index);
        omega_node.set(&omega_value, &omega_index);

        model.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);
        model.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);

        success *= (model.allocate() == 0);
        pmech_node.init(pmech0);
        success *= (model.verify() == 0);
        success *= (model.initialize() == 0);
        success *= (model.evaluateResidual() == 0);

        const ScalarT pmech_component =
            pmech0 * scalar(kConversionSystemBase / kConversionTrate);
        const ScalarT xflow0 = pmech_component + scalar(kDturb) * omega_value;
        const ScalarT vtemp0 = scalar(kAt) + scalar(kKt) * (scalar(kAt) - xflow0);

        success *= isEqual(model.y().getData()[index(Var::XVALVE)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::XFLOW)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::XTEMP)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VLOAD)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VTEMP)], vtemp0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VLV)], xflow0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::PMECH)], pmech0, scalar(kTolerance));

        checkZeroResidual(model, success);

        return success.report(__func__);
      }

      TestOutcome absoluteTolerance()
      {
        TestStatus success = true;

        Gov model(makeTestData());

        success                        *= (model.allocate() == 0);
        success                        *= (model.setAbsoluteTolerance(static_cast<RealT>(1.0e-7)) == 0);
        const auto& absolute_tolerance  = model.absoluteTolerance();
        success                        *= (absolute_tolerance.getSize() == static_cast<IdxT>(Var::MAXIMUM));

        const auto* tolerances = absolute_tolerance.getData();
        for (IdxT i = 0; i < absolute_tolerance.getSize(); ++i)
        {
          success *= isEqual(tolerances[i], scalar(1.0e-7), scalar(kTolerance));
        }

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
        ScalarT                                   pref_value{0.0};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      omega_index = static_cast<IdxT>(Var::MAXIMUM);
        IdxT                                      pref_index  = omega_index + 1;
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
        success *= isEqual(pref_node.read(),
                           prefForInitialPoint(pmech0, omega_value),
                           scalar(kTolerance));
        success *= (model.evaluateResidual() == 0);
        checkZeroResidual(model, success);

        pref_value += scalar(kPrefStep);
        success    *= (model.evaluateResidual() == 0);
        success    *= isEqual(model.getResidual().getData()[index(Var::VLOAD)],
                           scalar(kR * kPrefStep),
                           scalar(kTolerance));

        return success.report(__func__);
      }

      TestOutcome prefSignalBaseConversion()
      {
        TestStatus success = true;

        auto data                      = makeTestData();
        data.parameters[Params::Trate] = static_cast<RealT>(kConversionTrate);
        Gov model(data);
        model.setSystemBase(static_cast<RealT>(kSystemFrequency),
                            static_cast<RealT>(kConversionSystemBase * 1.0e6));

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        const ScalarT                             pmech0 = scalar(kConversionInitialPmech);
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   omega_value = scalar(kInitialOmega);
        ScalarT                                   pref_value{0.0};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      omega_index = static_cast<IdxT>(Var::MAXIMUM);
        IdxT                                      pref_index  = omega_index + 1;
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

        const ScalarT pmech_component =
            pmech0 * scalar(kConversionSystemBase / kConversionTrate);
        const ScalarT xflow0         = pmech_component + scalar(kDturb) * omega_value;
        const ScalarT pref_component = xflow0 + omega_value / scalar(kR);
        const ScalarT pref_system =
            pref_component * scalar(kConversionTrate / kConversionSystemBase);

        success *= isEqual(pref_node.read(), pref_system, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VLOAD)], xflow0, scalar(kTolerance));
        success *= (model.evaluateResidual() == 0);
        checkZeroResidual(model, success);

        pref_value += scalar(kPrefStep);
        success    *= (model.evaluateResidual() == 0);
        success    *= isEqual(model.getResidual().getData()[index(Var::VLOAD)],
                           scalar(kR * kConversionSystemBase / kConversionTrate * kPrefStep),
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
        IdxT                                      omega_index = static_cast<IdxT>(Var::MAXIMUM);
        IdxT                                      pref_index  = omega_index + 1;
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

        model.y().getData()[index(Var::XVALVE)] = xvalve;
        model.y().getData()[index(Var::XFLOW)]  = xflow;
        model.y().getData()[index(Var::XTEMP)]  = xtemp;
        model.y().getData()[index(Var::VLOAD)]  = vload;
        model.y().getData()[index(Var::VTEMP)]  = vtemp;
        model.y().getData()[index(Var::VLV)]    = vlv;
        model.y().getData()[index(Var::PMECH)]  = pmech;

        model.yp().getData()[index(Var::XVALVE)] = xvalve_dot;
        model.yp().getData()[index(Var::XFLOW)]  = xflow_dot;
        model.yp().getData()[index(Var::XTEMP)]  = xtemp_dot;

        model.y().setDataUpdated();
        model.yp().setDataUpdated();

        success *= (model.verify() == 0);
        success *= (model.evaluateResidual() == 0);

        const ScalarT              valve_target = vlv - xvalve;
        const ScalarT              selected_vlv = vload;
        const std::vector<ScalarT> expected     = {
            -xvalve_dot + valve_target / scalar(kT1),
            -xflow_dot + (-xflow + xvalve) / scalar(kT2),
            -xtemp_dot + (-xtemp + xflow) / scalar(kT3),
            -omega_value + scalar(kR) * (pref_value - vload),
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

        const RealT above_vmax = kVmax + kValveExcess;
        const RealT below_vmin = kVmin - kValveExcess;

        // Opening drive against the saturated upper limit is blocked.
        checkValveResidual(model, above_vmax, above_vmax + kValveDrive, 0.0, success);
        // Closing drive away from the saturated upper limit passes.
        checkValveResidual(model, above_vmax, above_vmax - kValveDrive, -kValveDrive / kT1, success);
        // Closing drive against the saturated lower limit is blocked.
        checkValveResidual(model, below_vmin, below_vmin - kValveDrive, 0.0, success);
        // Opening drive away from the saturated lower limit passes.
        checkValveResidual(model, below_vmin, below_vmin + kValveDrive, kValveDrive / kT1, success);

        return success.report(__func__);
      }

      TestOutcome responseModes()
      {
        TestStatus success = true;

        auto fixed_data                     = makeTestData();
        fixed_data.parameters[Params::mode] = static_cast<IdxT>(Mode::Fixed);
        fixed_data.parameters[Params::Vmin] = static_cast<RealT>(kResidualXvalve);
        fixed_data.parameters[Params::Vmax] = static_cast<RealT>(kResidualXvalve);

        Gov fixed_model(fixed_data);
        success *= (fixed_model.verify() == 0);
        success *= (fixed_model.allocate() == 0);
        setDynamicProbe(fixed_model);
        success *= (fixed_model.evaluateResidual() == 0);

        success *= isEqual(fixed_model.getResidual().getData()[index(Var::XVALVE)],
                           -scalar(kResidualXvalveDot),
                           scalar(kTolerance));
        success *= isEqual(fixed_model.getResidual().getData()[index(Var::XFLOW)],
                           -scalar(kResidualXflowDot),
                           scalar(kTolerance));
        success *= isEqual(fixed_model.getResidual().getData()[index(Var::XTEMP)],
                           -scalar(kResidualXtempDot),
                           scalar(kTolerance));

        auto down_only_data                     = makeTestData();
        down_only_data.parameters[Params::mode] = static_cast<IdxT>(Mode::DownOnly);

        Gov down_only_model(down_only_data);
        success *= (down_only_model.verify() == 0);
        success *= (down_only_model.allocate() == 0);
        setDynamicProbe(down_only_model);
        success *= (down_only_model.evaluateResidual() == 0);

        success *= isEqual(down_only_model.getResidual().getData()[index(Var::XVALVE)],
                           -scalar(kResidualXvalveDot)
                               + scalar(kResidualVlv - kResidualXvalve) / scalar(kT1),
                           scalar(kTolerance));
        success *= isEqual(down_only_model.getResidual().getData()[index(Var::XFLOW)],
                           -scalar(kResidualXflowDot)
                               + scalar(-kResidualXflow + kResidualXvalve) / scalar(kT2),
                           scalar(kTolerance));
        success *= isEqual(down_only_model.getResidual().getData()[index(Var::XTEMP)],
                           -scalar(kResidualXtempDot)
                               + scalar(-kResidualXtemp + kResidualXflow) / scalar(kT3),
                           scalar(kTolerance));

        return success.report(__func__);
      }

      TestOutcome narrowLimitWarning()
      {
        TestStatus success = true;

        using Log = Utilities::Logger;

        const auto         verbosity = Log::verbosity();
        std::ostringstream messages;
        Log::setOutput(messages);
        Log::setVerbosity(Log::WARNINGS);

        auto narrow_data                     = makeTestData();
        narrow_data.parameters[Params::Vmax] = static_cast<RealT>(0.01);
        Gov narrow_model(narrow_data);
        success *= (narrow_model.verify() == 0);
        success *= (messages.str().find("maximum anti-windup gate") != std::string::npos);

        messages.str("");
        messages.clear();
        Gov wide_model(makeTestData());
        success *= (wide_model.verify() == 0);
        success *= (messages.str().find("maximum anti-windup gate") == std::string::npos);

        messages.str("");
        messages.clear();
        auto fixed_data                     = narrow_data;
        fixed_data.parameters[Params::mode] = static_cast<IdxT>(Mode::Fixed);
        Gov fixed_model(fixed_data);
        success *= (fixed_model.verify() == 0);
        success *= (messages.str().find("maximum anti-windup gate") == std::string::npos);

        Log::setOutput(std::cout);
        Log::setVerbosity(verbosity);

        return success.report(__func__);
      }

      TestOutcome initializationValidation()
      {
        TestStatus success = true;

        auto valve_limited_data                   = makeTestData();
        // Keep the temperature margin open so only the valve limit is exercised.
        valve_limited_data.parameters[Params::At] = static_cast<RealT>(kVmax + 2.0 * kValveExcess);
        Gov valve_limited_model(valve_limited_data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> valve_limited_pmech_node;
        ScalarT                                   valve_limited_pmech_value{0.0};
        IdxT                                      valve_limited_pmech_index = INVALID_INDEX<IdxT>;
        valve_limited_pmech_node.set(&valve_limited_pmech_value, &valve_limited_pmech_index);

        valve_limited_model.getSignals().template assignSignalNode<Var::PMECH>(
            &valve_limited_pmech_node);

        success *= (valve_limited_model.allocate() == 0);
        valve_limited_pmech_node.init(scalar(kVmax + kValveExcess));
        // Over-rated dispatch warns but still initializes and evaluates.
        success *= (valve_limited_model.initialize() == 0);
        success *= (valve_limited_model.evaluateResidual() == 0);
        checkZeroResidual(valve_limited_model, success);

        auto temperature_limited_data                   = makeTestData();
        temperature_limited_data.parameters[Params::At] = static_cast<RealT>(0.0);
        Gov temperature_limited_model(temperature_limited_data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> temperature_limited_pmech_node;
        ScalarT                                   temperature_limited_pmech_value{0.0};
        IdxT                                      temperature_limited_pmech_index = INVALID_INDEX<IdxT>;
        temperature_limited_pmech_node.set(&temperature_limited_pmech_value,
                                           &temperature_limited_pmech_index);

        temperature_limited_model.getSignals().template assignSignalNode<Var::PMECH>(
            &temperature_limited_pmech_node);

        success *= (temperature_limited_model.allocate() == 0);
        temperature_limited_pmech_node.init(scalar(kInitialPmech));
        success *= (temperature_limited_model.initialize() != 0);

        return success.report(__func__);
      }

      TestOutcome smoothMinimumInitialization()
      {
        TestStatus success = true;

        const RealT at_value = kInitialPmech + kNearGateMargin / (1.0 + kKt);

        auto data                   = makeTestData();
        data.parameters[Params::At] = at_value;
        Gov model(data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        ScalarT                                   pmech_value{0.0};
        ScalarT                                   pref_value{0.0};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        IdxT                                      pref_index  = static_cast<IdxT>(Var::MAXIMUM);
        pmech_node.set(&pmech_value, &pmech_index);
        pref_node.set(&pref_value, &pref_index);

        model.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);
        model.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);

        success *= (model.allocate() == 0);
        pmech_node.init(scalar(kInitialPmech));
        success *= (model.initialize() == 0);
        success *= (model.evaluateResidual() == 0);

        const ScalarT xflow0 = scalar(kInitialPmech);
        const ScalarT at     = scalar(at_value);
        const ScalarT vtemp0 = at + scalar(kKt) * (at - xflow0);
        const RealT   margin = static_cast<RealT>(vtemp0 - xflow0);
        const ScalarT vload0 = vtemp0 - scalar(Math::iramp(margin));

        success *= (margin > ZERO<RealT>);
        success *= isEqual(margin, kNearGateMargin, kTolerance);
        success *= (vload0 > vtemp0);
        success *= isEqual(model.y().getData()[index(Var::VLOAD)], vload0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VTEMP)], vtemp0, scalar(kTolerance));
        success *= isEqual(model.y().getData()[index(Var::VLV)], xflow0, scalar(kTolerance));
        success *= isEqual(pref_node.read(), vload0, scalar(kTolerance));

        checkZeroResidual(model, success);

        return success.report(__func__);
      }

      TestOutcome smoothMinimumEqualityRejected()
      {
        TestStatus success = true;

        auto data                   = makeTestData();
        data.parameters[Params::At] = static_cast<RealT>(kInitialPmech);
        Gov model(data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech_node;
        ScalarT                                   pmech_value{0.0};
        IdxT                                      pmech_index = INVALID_INDEX<IdxT>;
        pmech_node.set(&pmech_value, &pmech_index);

        model.getSignals().template assignSignalNode<Var::PMECH>(&pmech_node);

        success *= (model.allocate() == 0);
        pmech_node.init(scalar(kInitialPmech));
        success *= (model.initialize() != 0);

        return success.report(__func__);
      }

      TestOutcome timeConstantMinimum()
      {
        TestStatus success = true;

        auto data                   = makeTestData();
        data.parameters[Params::T1] = static_cast<RealT>(0.0);
        data.parameters[Params::T2] = static_cast<RealT>(0.0);
        data.parameters[Params::T3] = static_cast<RealT>(0.0);

        Gov model(data);
        success *= (model.allocate() == 0);
        success *= (model.tagDifferentiable() == 0);
        success *= (model.tag()[index(Var::XVALVE)] == true);
        success *= (model.tag()[index(Var::XFLOW)] == true);
        success *= (model.tag()[index(Var::XTEMP)] == true);

        const ScalarT xvalve    = scalar(kResidualXvalve);
        const ScalarT xflow     = scalar(kResidualXflow);
        const ScalarT xtemp     = scalar(kResidualXtemp);
        const ScalarT vlv       = scalar(kResidualVlv);
        const ScalarT xflow_dot = scalar(kResidualXflowDot);
        const ScalarT probe_dot = scalar(kProbeDerivative);

        model.y().getData()[index(Var::XVALVE)] = xvalve;
        model.y().getData()[index(Var::XFLOW)]  = xflow;
        model.y().getData()[index(Var::XTEMP)]  = xtemp;
        model.y().getData()[index(Var::VLV)]    = vlv;

        model.yp().getData()[index(Var::XVALVE)] = probe_dot;
        model.yp().getData()[index(Var::XFLOW)]  = xflow_dot;
        model.yp().getData()[index(Var::XTEMP)]  = probe_dot;

        model.y().setDataUpdated();
        model.yp().setDataUpdated();

        success *= (model.evaluateResidual() == 0);
        success *= isEqual(model.getResidual().getData()[index(Var::XVALVE)],
                           -probe_dot + (vlv - xvalve) / scalar(kTimeConstantMinimum),
                           scalar(kTolerance));
        success *= isEqual(model.getResidual().getData()[index(Var::XFLOW)],
                           -xflow_dot + (-xflow + xvalve) / scalar(kTimeConstantMinimum),
                           scalar(kTolerance));
        success *= isEqual(model.getResidual().getData()[index(Var::XTEMP)],
                           -probe_dot + (-xtemp + xflow) / scalar(kTimeConstantMinimum),
                           scalar(kTolerance));

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        Gov default_model{Data{}};
        success *= (default_model.verify() == 0);

        auto zero_droop                  = makeTestData();
        zero_droop.parameters[Params::R] = ZERO<RealT>;
        Gov zero_droop_model(zero_droop);
        success *= (zero_droop_model.verify() > 0);

        auto negative_load_limit                   = makeTestData();
        negative_load_limit.parameters[Params::At] = static_cast<RealT>(-kAt);
        Gov negative_load_limit_model(negative_load_limit);
        success *= (negative_load_limit_model.verify() > 0);

        auto negative_temperature_gain                   = makeTestData();
        negative_temperature_gain.parameters[Params::Kt] = static_cast<RealT>(-kKt);
        Gov negative_temperature_gain_model(negative_temperature_gain);
        success *= (negative_temperature_gain_model.verify() > 0);

        auto negative_damping                      = makeTestData();
        negative_damping.parameters[Params::Dturb] = static_cast<RealT>(-kDturb);
        Gov negative_damping_model(negative_damping);
        success *= (negative_damping_model.verify() > 0);

        auto swapped_limits                     = makeTestData();
        swapped_limits.parameters[Params::Vmin] = static_cast<RealT>(kVmax);
        swapped_limits.parameters[Params::Vmax] = static_cast<RealT>(kVmin);
        Gov swapped_limits_model(swapped_limits);
        success *= (swapped_limits_model.verify() > 0);

        auto equal_normal_limits                     = makeTestData();
        equal_normal_limits.parameters[Params::Vmin] = static_cast<RealT>(kVmax);
        Gov equal_normal_limits_model(equal_normal_limits);
        success *= (equal_normal_limits_model.verify() > 0);

        auto invalid_mode                     = makeTestData();
        invalid_mode.parameters[Params::mode] = static_cast<IdxT>(3);
        Gov invalid_mode_model(invalid_mode);
        success *= (invalid_mode_model.verify() > 0);

        auto zero_trate                      = makeTestData();
        zero_trate.parameters[Params::Trate] = ZERO<RealT>;
        Gov zero_trate_model(zero_trate);
        success *= (zero_trate_model.verify() > 0);

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

      TestOutcome jsonParseAndSystemAssembly()
      {
        TestStatus success = true;

        std::istringstream input(R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "gas-turbine governor",
    "case_description": "GASTPTI parser and assembly test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    {
      "number": 1,
      "class": "infinite_bus",
      "name": "Bus 1",
      "init": { "Vr": 1.0, "Vi": 0.0 },
      "params": { "kv": 1.0 }
    }
  ],
  "signals": [
    { "signal_id": 10, "name": "Pmech" }
  ],
  "devices": [
    {
      "class": "Genrou",
      "ports": { "bus": 1, "pmech": 10 },
      "id": "GEN1",
      "params": {
        "p0": 0.3, "q0": 0.0, "H": 3.0, "D": 0.0, "Ra": 0.0,
        "Tdop": 7.0, "Tdopp": 0.04, "Tqop": 0.75, "Tqopp": 0.05,
        "Xd": 2.1, "Xdp": 0.2, "Xdpp": 0.18, "Xq": 0.5, "Xqp": 0.5,
        "Xqpp": 0.18, "Xl": 0.15, "S10": 0.0, "S12": 0.0,
        "mva": 100.0
      }
    },
    {
      "class": "GastPti",
      "ports": { "pmech": 10 },
      "id": "GAST1",
      "params": {
        "R": 0.05, "T1": 0.4, "T2": 0.5, "T3": 0.25,
        "At": 2.0, "Kt": 0.3, "Vmax": 0.3, "Vmin": 0.3,
        "Dturb": 0.1, "Trate": 100.0, "mode": 1
      },
      "mon": ["pmech", "xvalve"]
    }
  ]
}
)json");

        auto data             = PhasorDynamics::parseSystemModelData(input);
        success              *= (data.gastpti.size() == 1);
        const auto& gov_data  = data.gastpti[0];
        success              *= (gov_data.device_class == "GastPti");
        success              *= gov_data.buses.empty();
        success              *= gov_data.signal_inputs.empty();
        success              *= (gov_data.signal_outputs.at(Data::SignalOutputs::pmech)
                    == static_cast<IdxT>(10));
        success              *= (std::get_if<double>(&gov_data.parameters.at(Params::R))
                    != nullptr);
        success              *= (std::get_if<double>(&gov_data.parameters.at(Params::Trate))
                    != nullptr);
        success              *= (std::get_if<IdxT>(&gov_data.parameters.at(Params::mode))
                    != nullptr);
        success              *= (std::get<IdxT>(gov_data.parameters.at(Params::mode))
                    == static_cast<IdxT>(Mode::Fixed));
        success              *= (gov_data.monitored_variables.count(Mon::pmech) == 1);
        success              *= (gov_data.monitored_variables.count(Mon::xvalve) == 1);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);

        auto* governor  = dynamic_cast<Gov*>(system.getComponent(1));
        success        *= (governor != nullptr);
        if (governor != nullptr)
        {
          success *= isEqual(governor->y().getData()[index(Var::PMECH)],
                             scalar(0.3),
                             scalar(kSystemTolerance));
          const auto* pmech =
              governor->getSignals().template getSignalNode<Var::PMECH>();
          success *= (pmech != nullptr);
          if (pmech != nullptr)
          {
            success *= pmech->linked();
          }
        }

        success              *= (system.tagDifferentiable() == 0);
        success              *= (system.evaluateResidual() == 0);
        const auto& residual  = system.getResidual();
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          success *= isEqual(residual.getData()[i],
                             scalar(0.0),
                             scalar(kSystemTolerance));
        }
        success *= (system.evaluateJacobian() == 0);

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
      /// Component expectations computed in closed form match to machine precision.
      static constexpr RealT kTolerance       = std::numeric_limits<RealT>::epsilon();
      /// Smooth limiter and gate approximations only match their piecewise
      /// targets to the CommonMath smoothness width away from transitions.
      static constexpr RealT kSmoothTolerance = 1.0e-2;
      /// System assembly couples machine and governor initialization, which
      /// accumulates roundoff above machine precision.
      static constexpr RealT kSystemTolerance = 100.0 * std::numeric_limits<RealT>::epsilon();

      static constexpr RealT kR                   = 0.05;
      static constexpr RealT kT1                  = 0.4;
      static constexpr RealT kT2                  = 0.5;
      static constexpr RealT kT3                  = 0.25;
      static constexpr RealT kTimeConstantMinimum = 1.0e-3;
      static constexpr RealT kAt                  = 2.0;
      static constexpr RealT kKt                  = 0.3;
      static constexpr RealT kVmax                = 1.2;
      static constexpr RealT kVmin                = 0.0;
      static constexpr RealT kDturb               = 0.1;
      static constexpr RealT kTrate               = 100.0;

      static constexpr RealT kInitialPmech   = 0.75;
      static constexpr RealT kInitialOmega   = 0.02;
      static constexpr RealT kNearGateMargin = 1.0e-4;
      static constexpr RealT kPrefStep       = 0.1;

      /// Distance outside the valve limits used to saturate the anti-windup.
      static constexpr RealT kValveExcess     = 1.0;
      /// Magnitude of the valve drive applied against a saturated limit.
      static constexpr RealT kValveDrive      = 1.0;
      /// Derivative value verifying floored lags keep their derivative terms.
      static constexpr RealT kProbeDerivative = 5.0;

      static constexpr RealT kSystemFrequency        = 60.0;
      static constexpr RealT kConversionTrate        = 50.0;
      static constexpr RealT kConversionSystemBase   = 100.0;
      static constexpr RealT kConversionInitialPmech = 0.40;

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

      template <typename value_type>
      static value_type value(RealT value)
      {
        return value_type{value};
      }

      static size_t index(Var variable)
      {
        return static_cast<size_t>(variable);
      }

      template <typename value_type>
      static value_type prefForInitialPoint(const value_type& pmech, const value_type& omega)
      {
        return pmech + value<value_type>(kDturb) * omega + omega / value<value_type>(kR);
      }

      void checkResidual(const Gov&                  model,
                         const std::vector<ScalarT>& expected,
                         TestStatus&                 success) const
      {
        const auto& residual = model.getResidual();
        for (size_t i = 0; i < expected.size(); ++i)
        {
          if (!isEqual(residual.getData()[i], expected[i], scalar(kTolerance)))
          {
            std::cout << "Unexpected GASTPTI residual at index " << i << ": "
                      << residual.getData()[i] << " != " << expected[i] << "\n";
            success = false;
          }
        }
      }

      void checkValveResidual(Gov&        model,
                              RealT       xvalve,
                              RealT       vlv,
                              RealT       expected,
                              TestStatus& success) const
      {
        model.y().setToZero();
        model.yp().setToZero();
        model.y().getData()[index(Var::XVALVE)] = scalar(xvalve);
        model.y().getData()[index(Var::VLV)]    = scalar(vlv);
        model.y().setDataUpdated();
        success *= (model.evaluateResidual() == 0);
        success *= isEqual(model.getResidual().getData()[index(Var::XVALVE)],
                           scalar(expected),
                           scalar(kSmoothTolerance / kT1));
      }

      void setDynamicProbe(Gov& model) const
      {
        model.y().setToZero();
        model.yp().setToZero();

        model.y().getData()[index(Var::XVALVE)] = scalar(kResidualXvalve);
        model.y().getData()[index(Var::XFLOW)]  = scalar(kResidualXflow);
        model.y().getData()[index(Var::XTEMP)]  = scalar(kResidualXtemp);
        model.y().getData()[index(Var::VLV)]    = scalar(kResidualVlv);

        model.yp().getData()[index(Var::XVALVE)] = scalar(kResidualXvalveDot);
        model.yp().getData()[index(Var::XFLOW)]  = scalar(kResidualXflowDot);
        model.yp().getData()[index(Var::XTEMP)]  = scalar(kResidualXtempDot);

        model.y().setDataUpdated();
        model.yp().setDataUpdated();
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
        data.monitored_variables.insert(Mon::xvalve);

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
        data.parameters[Params::mode]  = static_cast<IdxT>(Mode::Normal);

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

        model.y().getData()[index(Var::XVALVE)] = ADScalarT{kResidualXvalve};
        model.y().getData()[index(Var::XFLOW)]  = ADScalarT{kResidualXflow};
        model.y().getData()[index(Var::XTEMP)]  = ADScalarT{kResidualXtemp};
        model.y().getData()[index(Var::VLOAD)]  = ADScalarT{kResidualVload};
        model.y().getData()[index(Var::VTEMP)]  = ADScalarT{kResidualVtemp};
        model.y().getData()[index(Var::VLV)]    = ADScalarT{kResidualVlv};
        model.y().getData()[index(Var::PMECH)]  = ADScalarT{kResidualPmech};

        model.yp().getData()[index(Var::XVALVE)] = ADScalarT{kResidualXvalveDot};
        model.yp().getData()[index(Var::XFLOW)]  = ADScalarT{kResidualXflowDot};
        model.yp().getData()[index(Var::XTEMP)]  = ADScalarT{kResidualXtempDot};

        for (size_t i = 0; i < model.y().getSize(); ++i)
        {
          model.y().getData()[i].setVariableNumber(i);
          model.yp().getData()[i].setVariableNumber(i);
        }
        model.y().setDataUpdated();
        model.yp().setDataUpdated();
        omega_value.setVariableNumber(static_cast<size_t>(omega_index));
        pref_value.setVariableNumber(static_cast<size_t>(pref_index));

        model.evaluateResidual();

        std::vector<DependencyMap> dependencies(model.getResidual().getSize());
        for (size_t i = 0; i < dependencies.size(); ++i)
        {
          dependencies[i] = model.getResidual().getData()[i].getDependencies();
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

        model.y().getData()[index(Var::XVALVE)] = scalar(kResidualXvalve);
        model.y().getData()[index(Var::XFLOW)]  = scalar(kResidualXflow);
        model.y().getData()[index(Var::XTEMP)]  = scalar(kResidualXtemp);
        model.y().getData()[index(Var::VLOAD)]  = scalar(kResidualVload);
        model.y().getData()[index(Var::VTEMP)]  = scalar(kResidualVtemp);
        model.y().getData()[index(Var::VLV)]    = scalar(kResidualVlv);
        model.y().getData()[index(Var::PMECH)]  = scalar(kResidualPmech);

        model.yp().getData()[index(Var::XVALVE)] = scalar(kResidualXvalveDot);
        model.yp().getData()[index(Var::XFLOW)]  = scalar(kResidualXflowDot);
        model.yp().getData()[index(Var::XTEMP)]  = scalar(kResidualXtempDot);

        model.y().setDataUpdated();
        model.yp().setDataUpdated();

        model.evaluateResidual();
        model.evaluateJacobian();

        model.constructCsr();
        auto* model_jacobian = model.getCsrJacobian();

        return MapFromCsr(model_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
