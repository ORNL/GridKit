#pragma once

#include <array>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    class ControllerRepcaTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ControllerRepcaTests()  = default;
      ~ControllerRepcaTests() = default;

      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

      /// Validate construction, defaults, parameters, signals, and time floors.
      TestOutcome validation()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors and warnings from the invalid cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Controller::Repca<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(Utilities::enum_size<Vars>()));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() > 0);

        Fixture<ScalarT> configured(makeData());
        configured.attachAllInputs();
        success *= (configured.repca.size() == static_cast<IdxT>(Utilities::enum_size<Vars>()));
        success *= (configured.repca.getMonitor() != nullptr);
        success *= (configured.repca.verify() == 0);

        Fixture<ScalarT> documented_defaults(makeMinimalData());
        documented_defaults.attachAllInputs();
        success *= (documented_defaults.repca.verify() == 0);
        success *= defaultsMatchDocumentedValues();

        auto integer_numeric                    = makeData();
        integer_numeric.parameters[Params::mva] = static_cast<IdxT>(100);
        Fixture<ScalarT> integer_parameter(integer_numeric);
        integer_parameter.attachAllInputs();
        success *= (integer_parameter.repca.verify() == 0);

        PhasorDynamics::Controller::Repca<ScalarT, IdxT> missing_signals(&bus, makeData());
        success *= (missing_signals.verify() > 0);

        success *= invalidParameterCase(Params::mva, 0.0);
        success *= invalidParameterCase(Params::Tfv, -0.1);
        success *= invalidParameterCase(Params::dbdlow, 0.1);
        success *= invalidParameterCase(Params::dbdupper, -0.1);
        success *= invalidParameterCase(Params::emin, 0.1);
        success *= invalidParameterCase(Params::emax, -0.1);
        success *= invalidParameterCase(Params::Qmin, 1.1);
        success *= invalidParameterCase(Params::fdbd1, 0.1);
        success *= invalidParameterCase(Params::fdbd2, -0.1);
        success *= invalidParameterCase(Params::Ddn, -0.1);
        success *= invalidParameterCase(Params::Dup, -0.1);
        success *= invalidParameterCase(Params::femin, 0.1);
        success *= invalidParameterCase(Params::femax, -0.1);
        success *= invalidParameterCase(Params::Pmin, 2.1);
        success *= invalidParameterCase(Params::mva, true);

        success *= invalidParameterCase(Params::Tfltr, -0.2);
        success *= invalidParameterCase(Params::Tft, -0.1);
        success *= invalidParameterCase(Params::Tp, -0.3);
        success *= invalidParameterCase(Params::Tlag, -0.4);

        const RealT                  nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT                  infinity = std::numeric_limits<RealT>::infinity();
        const std::array<Params, 28> real_parameters{{
            Params::mva,
            Params::Tfltr,
            Params::Vfrz,
            Params::Rc,
            Params::Xc,
            Params::Kc,
            Params::dbdlow,
            Params::dbdupper,
            Params::emax,
            Params::emin,
            Params::Kp,
            Params::Ki,
            Params::Qmax,
            Params::Qmin,
            Params::Tft,
            Params::Tfv,
            Params::Tp,
            Params::fdbd1,
            Params::fdbd2,
            Params::Ddn,
            Params::Dup,
            Params::femax,
            Params::femin,
            Params::Kpg,
            Params::Kig,
            Params::Pmax,
            Params::Pmin,
            Params::Tlag,
        }};
        for (const Params parameter : real_parameters)
        {
          success *= invalidParameterCase(parameter, nan);
          success *= invalidParameterCase(parameter, infinity);
          success *= invalidParameterCase(parameter, -infinity);
        }
        success *= invalidParameterCase(Params::mva, std::numeric_limits<RealT>::max());

        {
          Fixture<ScalarT> nonfinite_system_base(makeData(), 1.0, 0.0, infinity);
          nonfinite_system_base.attachAllInputs();
          success *= (nonfinite_system_base.repca.verify() > 0);
        }
        {
          auto tiny_base_data                    = makeData();
          tiny_base_data.parameters[Params::mva] = std::numeric_limits<RealT>::min();
          Fixture<ScalarT> overflowing_base_ratio(tiny_base_data,
                                                  1.0,
                                                  0.0,
                                                  std::numeric_limits<RealT>::max());
          overflowing_base_ratio.attachAllInputs();
          success *= (overflowing_base_ratio.repca.verify() > 0);
        }

        const std::array<Params, 3> flag_parameters{{
            Params::VcompFlag,
            Params::RefFlag,
            Params::Freqflag,
        }};
        const std::array<bool, 2>   valid_flag_values{{false, true}};
        const std::array<IdxT, 3>   invalid_integral_flag_values{{
            static_cast<IdxT>(0),
            static_cast<IdxT>(1),
            static_cast<IdxT>(2),
        }};
        const std::array<RealT, 5>  invalid_real_flag_values{{
            static_cast<RealT>(0.0),
            static_cast<RealT>(0.5),
            static_cast<RealT>(1.0),
            nan,
            infinity,
        }};
        for (const Params flag : flag_parameters)
        {
          for (const bool value : valid_flag_values)
          {
            auto data             = makeData();
            data.parameters[flag] = value;
            Fixture<ScalarT> model(data);
            model.attachAllInputs();
            success *= (model.repca.verify() == 0);
          }

          for (const IdxT value : invalid_integral_flag_values)
          {
            success *= invalidParameterCase(flag, value);
          }

          for (const RealT value : invalid_real_flag_values)
          {
            success *= invalidParameterCase(flag, value);
          }
        }

        PhasorDynamics::Controller::Repca<ScalarT, IdxT> busless(nullptr, makeData());
        success *= (busless.verify() > 0);

        success *= unlinkedSignalRejected<Ext::ir>();
        success *= unlinkedSignalRejected<Ext::ii>();
        success *= unlinkedSignalRejected<Ext::p>();
        success *= unlinkedSignalRejected<Ext::q>();
        success *= unlinkedSignalRejected<Ext::freq>();
        success *= unlinkedSignalRejected<Ext::vref>();
        success *= unlinkedSignalRejected<Ext::pref>();
        success *= unlinkedSignalRejected<Ext::qref>();
        success *= unlinkedSignalRejected<Ext::freqref>();

        auto floor_data                      = makeInitializationData();
        floor_data.parameters[Params::Tfltr] = 0.0;
        floor_data.parameters[Params::Tfv]   = 0.0;
        floor_data.parameters[Params::Tp]    = 0.0;
        floor_data.parameters[Params::Tlag]  = 0.0;

        Fixture<ScalarT> floored(floor_data);
        floored.attachAllInputs();
        setInitializationInputs(floored);
        success *= floored.initialize(0.25, 0.45);
        success *= (floored.repca.evaluateResidual() == 0);
        success *= allResidualsWithinInitTolerance(floored.repca);

        auto* y                = floored.repca.y().getData();
        y[index(Vars::VMEAS)] -= 0.001;
        y[index(Vars::QMEAS)] -= 0.002;
        y[index(Vars::PMEAS)] -= 0.003;
        y[index(Vars::PREF)]  -= 0.004;
        floored.repca.y().setDataUpdated();
        success *= (floored.repca.evaluateResidual() == 0);
        // The PMEAS step also raises the active error, which the PI output
        // passes to the command lag: (0.004 + Kpg * 0.003) / 0.001.
        const std::array<VariableValue, 4> floored_residuals{{
            {Vars::VMEAS, 1.0},
            {Vars::QMEAS, 2.0},
            {Vars::PMEAS, 3.0},
            {Vars::PREF, 9.1},
        }};
        success *= residualsMatch(floored.repca,
                                  floored_residuals,
                                  "floored time constants");

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      /// Check initialization state, signal publication, monitors, and flag modes.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeInitializationData(), 0.8, 0.6);
        fixture.attachAllInputs(99.0);
        setInitializationInputs(fixture);
        success *= fixture.initialize(0.25, 0.45);
        success *= (fixture.repca.tagDifferentiable() == 0);
        success *= (fixture.repca.evaluateResidual() == 0);

        const std::array<VariableValue, Utilities::enum_size<Vars>()> initial_state{{
            {Vars::VMEAS, 0.984002032518226},
            {Vars::QMEAS, 0.2},
            {Vars::XQPI, 0.5},
            {Vars::XQLAG, 0.5},
            {Vars::PMEAS, 0.8},
            {Vars::XPPI, 0.9},
            {Vars::PREF, 0.9},
            {Vars::QEXT, 0.25},
            {Vars::PEXT, 0.45},
        }};
        success *= stateMatches(fixture.repca, initial_state, "initialization");

        success *= scalarPreserved(fixture.input(Ext::ir), 0.2, "preserved ir");
        success *= scalarPreserved(fixture.input(Ext::ii), -0.1, "preserved ii");
        success *= scalarPreserved(fixture.input(Ext::p), 0.4, "preserved p");
        success *= scalarPreserved(fixture.input(Ext::q), 0.1, "preserved q");
        success *= scalarPreserved(fixture.input(Ext::freq), 0.99, "preserved freq");
        success *= scalarPreserved(fixture.qext(), 0.25, "preserved qext");
        success *= scalarPreserved(fixture.pext(), 0.45, "preserved pext");
        success *= scalarMatches(fixture.input(Ext::vref), 0.984002032518226, "published vref");
        success *= scalarMatches(fixture.input(Ext::pref), 0.4, "published pref");
        success *= scalarMatches(fixture.input(Ext::qref), 0.1, "published qref");
        success *= scalarMatches(fixture.input(Ext::freqref), 0.99, "published freqref");

        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          const bool expected = row <= index(Vars::PREF);
          if (fixture.repca.tag()[row] != expected)
          {
            std::cout << "REPCA differentiability tag " << row << " mismatch\n";
            success = false;
          }
        }

        constexpr RealT absolute_tolerance = 2.5e-7;

        success *= (fixture.repca.setAbsoluteTolerance(absolute_tolerance) == 0);

        const auto* tolerances = fixture.repca.absoluteTolerance().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          success *= valueUnchanged(tolerances[row],
                                    absolute_tolerance,
                                    "absolute tolerance",
                                    row);
        }

        const std::array<RealT, 5> initial_monitors{{
            0.25,
            0.45,
            0.984002032518226,
            0.2,
            0.8,
        }};
        success *= monitorMatches(fixture.repca, initial_monitors, "initialization");

        success *= allResidualsWithinInitTolerance(fixture.repca);

        {
          auto exact_data                    = makeInitializationData();
          exact_data.parameters[Params::mva] = 73.0;

          // This non-binary base ratio exposes any component-base round trip.
          Fixture<ScalarT> exact_commands(exact_data, 0.8, 0.6);
          exact_commands.attachAllInputs();
          setInitializationInputs(exact_commands);
          success *= exact_commands.initialize(0.25, 0.45);
          success *= scalarPreserved(exact_commands.qext(), 0.25, "qext signal");
          success *= scalarPreserved(exact_commands.pext(), 0.45, "pext signal");
          success *= scalarPreserved(exact_commands.input(Ext::qref), 0.1, "qref signal");
          success *= (exact_commands.repca.evaluateResidual() == 0);
          success *= allResidualsWithinInitTolerance(exact_commands.repca);
        }

        Fixture<ScalarT> fallback(makeInitializationData(), 0.8, 0.6);
        fallback.attachAllInputs(0.0, false);
        setInitializationInputs(fallback);
        const auto previous_verbosity = Log::verbosity();
        // Suppress the expected missing-frequency warning for this fallback case.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        success *= fallback.initialize(0.25, 0.45);
        Log::setVerbosity(previous_verbosity);
        success *= scalarMatches(fallback.input(Ext::freqref), 1.0, "default frequency");
        success *= (fallback.repca.evaluateResidual() == 0);
        success *= allResidualsWithinInitTolerance(fallback.repca);

        Fixture<ScalarT> outputless(makeInitializationData(),
                                    0.8,
                                    0.6,
                                    100.0e6,
                                    false);
        outputless.attachAllInputs();
        setInitializationInputs(outputless);
        success *= outputless.initialize(0.25, 0.45);
        success *= (outputless.repca.evaluateResidual() == 0);
        const std::array<VariableValue, 2> outputless_state{{
            {Vars::QEXT, 0.25},
            {Vars::PEXT, 0.45},
        }};
        success *= stateMatches(outputless.repca,
                                outputless_state,
                                "unassigned command outputs");
        success *= monitorMatches(outputless.repca,
                                  initial_monitors,
                                  "unassigned command outputs");
        success *= allResidualsWithinInitTolerance(outputless.repca);

        struct FlagCase
        {
          const char* label;
          bool        voltage_compensation;
          bool        voltage_reference;
          bool        frequency_control;
          RealT       voltage;
          RealT       pref;
          RealT       pext;
        };

        const std::array<FlagCase, 8> flag_cases{{
            {"droop/reactive-reference/disabled-frequency", false, false, false, 1.08, 0.8, 0.0},
            {"droop/reactive-reference/enabled-frequency", false, false, true, 1.08, 0.9, 0.45},
            {"droop/voltage-reference/disabled-frequency", false, true, false, 1.08, 0.8, 0.0},
            {"droop/voltage-reference/enabled-frequency", false, true, true, 1.08, 0.9, 0.45},
            {"line-drop/reactive-reference/disabled-frequency", true, false, false, 0.984002032518226, 0.8, 0.0},
            {"line-drop/reactive-reference/enabled-frequency", true, false, true, 0.984002032518226, 0.9, 0.45},
            {"line-drop/voltage-reference/disabled-frequency", true, true, false, 0.984002032518226, 0.8, 0.0},
            {"line-drop/voltage-reference/enabled-frequency", true, true, true, 0.984002032518226, 0.9, 0.45},
        }};
        for (const auto& test_case : flag_cases)
        {
          auto data                          = makeInitializationData();
          data.parameters[Params::VcompFlag] = test_case.voltage_compensation;
          data.parameters[Params::RefFlag]   = test_case.voltage_reference;
          data.parameters[Params::Freqflag]  = test_case.frequency_control;

          Fixture<ScalarT> scenario(data, 0.8, 0.6);
          scenario.attachAllInputs(99.0);
          setInitializationInputs(scenario);
          success *= scenario.initialize(0.25, 0.45);
          success *= (scenario.repca.evaluateResidual() == 0);

          const std::array<VariableValue, 3> expected_state{{
              {Vars::VMEAS, test_case.voltage},
              {Vars::PREF, test_case.pref},
              {Vars::PEXT, test_case.pext},
          }};
          success *= stateMatches(scenario.repca,
                                  expected_state,
                                  test_case.label);
          success *= scalarMatches(scenario.qext(), 0.25, test_case.label);
          success *= scalarMatches(scenario.pext(), test_case.pext, test_case.label);
          success *= allResidualsWithinInitTolerance(scenario.repca);
        }

        return success.report(__func__);
      }

      /// Check initialization domains, adjusted limits, and atomicity.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors and limit-adjustment warnings from the cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        const auto data = makeInitializationData();

        struct RejectionCase
        {
          const char* label;
          RealT       qext;
          RealT       pext;
        };

        const RealT                        nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT                        infinity = std::numeric_limits<RealT>::infinity();
        const std::array<RejectionCase, 4> rejection_cases{{
            {"nonfinite qext", infinity, 0.45},
            {"nonfinite pext", 0.25, infinity},
            {"nan qext", nan, 0.45},
            {"nan pext", 0.25, nan},
        }};
        for (const auto& test_case : rejection_cases)
        {
          success *= initializationRejectedAtomically(data,
                                                      test_case.qext,
                                                      test_case.pext,
                                                      test_case.label);
        }
        const std::array<Ext, 5>   required_ports{{
            Ext::ir,
            Ext::ii,
            Ext::p,
            Ext::q,
            Ext::freq,
        }};
        const std::array<RealT, 3> nonfinite_values{{infinity, -infinity, nan}};
        for (const Ext port : required_ports)
        {
          for (const RealT value : nonfinite_values)
          {
            success *= initializationRejectedAtomically(data,
                                                        0.25,
                                                        0.45,
                                                        "nonfinite required signal",
                                                        NonfiniteTarget::INPUT,
                                                        0.8,
                                                        0.6,
                                                        port,
                                                        value);
          }
        }
        success *= initializationRejectedAtomically(data,
                                                    0.25,
                                                    0.45,
                                                    "nonfinite bus voltage",
                                                    NonfiniteTarget::BUS_VOLTAGE);

        auto collapsed_data = data;

        collapsed_data.parameters[Params::Qmin] = 0.5;
        collapsed_data.parameters[Params::Qmax] = 0.5;
        collapsed_data.parameters[Params::Pmin] = 0.9;
        collapsed_data.parameters[Params::Pmax] = 0.9;

        auto reactive_aw_data                         = makeInitializationData();
        reactive_aw_data.parameters[Params::dbdupper] = 0.03;
        reactive_aw_data.parameters[Params::emin]     = 0.0;

        // The published reference reproduces the selected reactive error,
        // whose deadbanded value takes the -0.1 offset branch of the limiter
        // inverse below emin = 0.
        const std::array<bool, 2>                    voltage_reference_values{{false, true}};
        const std::array<std::pair<RealT, RealT>, 2> voltage_cases{{
            {0.8, 0.6},
            {0.2, 0.0},
        }};
        for (const bool voltage_reference : voltage_reference_values)
        {
          auto data                        = reactive_aw_data;
          data.parameters[Params::RefFlag] = voltage_reference;
          for (const auto& voltage : voltage_cases)
          {
            Fixture<ScalarT> asymmetric_reactive(data,
                                                 voltage.first,
                                                 voltage.second);
            asymmetric_reactive.attachAllInputs();
            setInitializationInputs(asymmetric_reactive);
            success                    *= asymmetric_reactive.initialize(0.25, 0.45);
            success                    *= (asymmetric_reactive.repca.evaluateResidual() == 0);
            const RealT reactive_error  = reactiveError(asymmetric_reactive, voltage_reference);
            success                    *= scalarMatches(Math::deadband2(reactive_error, -0.02, 0.03),
                                     -0.1,
                                     "asymmetric reactive initialization");
            success                    *= allResidualsWithinInitTolerance(asymmetric_reactive.repca);
          }
        }

        auto frozen_data                     = reactive_aw_data;
        frozen_data.parameters[Params::Vfrz] = 0.9;

        {
          Fixture<ScalarT> frozen_reactive_rate(frozen_data, 0.05, 0.0);
          frozen_reactive_rate.attachAllInputs();
          setInitializationInputs(frozen_reactive_rate);
          success *= frozen_reactive_rate.initialize(0.25, 0.45);
          success *= (frozen_reactive_rate.repca.evaluateResidual() == 0);
          success *= allResidualsWithinInitTolerance(frozen_reactive_rate.repca);
        }

        {
          auto active_aw_data                      = makeInitializationData();
          active_aw_data.parameters[Params::fdbd2] = 0.025;
          active_aw_data.parameters[Params::femin] = 0.0;
          Fixture<ScalarT> asymmetric_active(active_aw_data, 0.8, 0.6);
          asymmetric_active.attachAllInputs();
          setInitializationInputs(asymmetric_active);
          success *= asymmetric_active.initialize(0.25, 0.45);
          success *= (asymmetric_active.repca.evaluateResidual() == 0);
          // The published plant reference carries the -0.1 lower-boundary error.
          success *= scalarMatches(asymmetric_active.input(Ext::freqref),
                                   0.995,
                                   "asymmetric frequency reference");
          success *= scalarMatches(asymmetric_active.input(Ext::pref),
                                   0.35,
                                   "lower-bound plant reference");
          success *= allResidualsWithinInitTolerance(asymmetric_active.repca);
        }

        {
          auto active_aw_data                      = makeInitializationData();
          active_aw_data.parameters[Params::fdbd1] = -0.025;
          active_aw_data.parameters[Params::femax] = 0.0;
          Fixture<ScalarT> asymmetric_active(active_aw_data, 0.8, 0.6);
          asymmetric_active.attachAllInputs();
          setInitializationInputs(asymmetric_active);
          success *= asymmetric_active.initialize(0.25, 0.45);
          success *= (asymmetric_active.repca.evaluateResidual() == 0);
          // The published plant reference carries the 0.1 upper-boundary error.
          success *= scalarMatches(asymmetric_active.input(Ext::freqref),
                                   0.985,
                                   "asymmetric frequency reference");
          success *= scalarMatches(asymmetric_active.input(Ext::pref),
                                   0.45,
                                   "upper-bound plant reference");
          success *= allResidualsWithinInitTolerance(asymmetric_active.repca);
        }

        auto overflow_data                    = makeInitializationData();
        overflow_data.parameters[Params::Rc]  = std::numeric_limits<RealT>::max();
        success                              *= initializationRejectedAtomically(overflow_data,
                                                    0.25,
                                                    0.45,
                                                    "nonfinite derived initialization candidate");

        // An invalid configuration is rejected before any state is written.
        {
          auto invalid_data                    = data;
          invalid_data.parameters[Params::Tfv] = -0.1;
          Fixture<ScalarT> invalid_fixture(invalid_data);
          invalid_fixture.attachAllInputs();
          setInitializationInputs(invalid_fixture);
          success *= (invalid_fixture.repca.allocate() == 0);
          poisonState(invalid_fixture, 0.25, 0.45);
          const auto invalid_y  = copyVector(invalid_fixture.repca.y());
          const auto invalid_yp = copyVector(invalid_fixture.repca.yp());
          if (invalid_fixture.repca.initialize() == 0)
          {
            std::cout << "Expected REPCA initialization rejection: invalid configuration\n";
            success = false;
          }
          success *= vectorUnchanged(invalid_fixture.repca.y(), invalid_y, "state");
          success *= vectorUnchanged(invalid_fixture.repca.yp(), invalid_yp, "derivative");
        }

        // A command exactly on a limit is reconstructed through the offset
        // branch of the limiter inverse, which leaves a smoothing-scaled
        // residual, so only the reconstructed state is checked.
        {
          Fixture<ScalarT> qmax_pmin_boundary(data, 0.8, 0.6);
          qmax_pmin_boundary.attachAllInputs();
          setInitializationInputs(qmax_pmin_boundary);
          success *= qmax_pmin_boundary.initialize(0.75, 0.0);
          success *= (qmax_pmin_boundary.repca.evaluateResidual() == 0);
          const std::array<VariableValue, 6> qmax_pmin_state{{
              {Vars::XQLAG, 1.5},
              {Vars::XQPI, 1.6},
              {Vars::QEXT, 0.75},
              {Vars::PREF, 0.0},
              {Vars::XPPI, -0.1},
              {Vars::PEXT, 0.0},
          }};
          success *= stateMatches(qmax_pmin_boundary.repca,
                                  qmax_pmin_state,
                                  "Qmax/Pmin command boundary");
        }

        {
          Fixture<ScalarT> qmin_pmax_boundary(data, 0.8, 0.6);
          qmin_pmax_boundary.attachAllInputs();
          setInitializationInputs(qmin_pmax_boundary);
          success *= qmin_pmax_boundary.initialize(-0.4, 1.0);
          success *= (qmin_pmax_boundary.repca.evaluateResidual() == 0);
          const std::array<VariableValue, 6> qmin_pmax_state{{
              {Vars::XQLAG, -0.8},
              {Vars::XQPI, -0.9},
              {Vars::QEXT, -0.4},
              {Vars::PREF, 2.0},
              {Vars::XPPI, 2.1},
              {Vars::PEXT, 1.0},
          }};
          success *= stateMatches(qmin_pmax_boundary.repca,
                                  qmin_pmax_state,
                                  "Qmin/Pmax command boundary");
        }

        {
          Fixture<ScalarT> collapsed_limits(collapsed_data, 0.8, 0.6);
          collapsed_limits.attachAllInputs();
          setInitializationInputs(collapsed_limits);
          success *= collapsed_limits.initialize(0.25, 0.45);
          success *= (collapsed_limits.repca.evaluateResidual() == 0);
          const std::array<VariableValue, 6> collapsed_state{{
              {Vars::XQPI, 0.5},
              {Vars::XQLAG, 0.5},
              {Vars::QEXT, 0.25},
              {Vars::XPPI, 0.9},
              {Vars::PREF, 0.9},
              {Vars::PEXT, 0.45},
          }};
          success *= stateMatches(collapsed_limits.repca,
                                  collapsed_state,
                                  "collapsed Q/P limits");
          success *= allResidualsWithinInitTolerance(collapsed_limits.repca);
        }

        struct LimitCase
        {
          const char* label;
          RealT       qext;
          RealT       pext;
          Vars        output;
          RealT       expected;
        };

        // At rest the reactive lag holds the PI output and the command lag
        // holds the active PI output.
        const std::array<LimitCase, 4> limit_cases{{
            {"adjusted Qmin", -0.5, 0.45, Vars::XQLAG, -1.0},
            {"adjusted Qmax", 0.9, 0.45, Vars::XQLAG, 1.8},
            {"adjusted Pmin", 0.25, -0.1, Vars::PREF, -0.2},
            {"adjusted Pmax", 0.25, 1.1, Vars::PREF, 2.2},
        }};

        for (const auto& test_case : limit_cases)
        {
          Fixture<ScalarT> adjusted(data, 0.8, 0.6);
          adjusted.attachAllInputs();
          setInitializationInputs(adjusted);
          success *= adjusted.initialize(test_case.qext, test_case.pext);
          success *= stateMatches(adjusted.repca,
                                  {{test_case.output, test_case.expected}},
                                  test_case.label);
          success *= (adjusted.repca.evaluateResidual() == 0);
          success *= allResidualsWithinInitTolerance(adjusted.repca);
        }

        {
          auto disabled_data                         = data;
          disabled_data.parameters[Params::Freqflag] = false;
          disabled_data.parameters[Params::Pmax]     = 0.5;

          Fixture<ScalarT> disabled_frequency(disabled_data, 0.8, 0.6);
          disabled_frequency.attachAllInputs();
          setInitializationInputs(disabled_frequency);
          success *= disabled_frequency.initialize(0.25, 0.45);
          success *= stateMatches(disabled_frequency.repca,
                                  {{Vars::PREF, 0.8},
                                   {Vars::PEXT, 0.0}},
                                  "measured-power limit");
          success *= (disabled_frequency.repca.evaluateResidual() == 0);
          success *= allResidualsWithinInitTolerance(disabled_frequency.repca);

          // A 0.05 system-base reference step is a 0.1 active error; the PI
          // state then places its output at 0.65 inside the 0.8 limit.
          disabled_frequency.input(Ext::pref) += 0.05;
          setState(disabled_frequency.repca, {{Vars::XPPI, 0.48}});
          setDerivative(disabled_frequency.repca, {{Vars::XPPI, 0.0}});
          success *= (disabled_frequency.repca.evaluateResidual() == 0);
          success *= residualsMatch(disabled_frequency.repca,
                                    {{Vars::XPPI, 0.18}},
                                    "measured-power limit");
        }

        {
          Fixture<ScalarT> adjusted(data, 0.8, 0.6);
          adjusted.attachAllInputs();
          setInitializationInputs(adjusted);
          success                    *= adjusted.initialize(1.0, 1.25);
          const RealT published_vref  = adjusted.input(Ext::vref);
          const RealT published_pref  = adjusted.input(Ext::pref);

          // Reference steps give 0.2 reactive and active errors; the PI states
          // place the outputs at 1.75 and 2.25, inside the adjusted limits.
          adjusted.input(Ext::vref) = published_vref + 0.22;
          adjusted.input(Ext::pref) = published_pref + 0.1;
          setState(adjusted.repca, {{Vars::XQPI, 1.35}, {Vars::XPPI, 1.91}});
          setDerivative(adjusted.repca, {{Vars::XQPI, 0.0}, {Vars::XPPI, 0.0}});
          success *= (adjusted.repca.evaluateResidual() == 0);
          success *= residualsMatch(adjusted.repca,
                                    {{Vars::XQPI, 0.6}, {Vars::XPPI, 0.36}},
                                    "adjusted antiwindup limits");

          // With no error, the lags read the PI outputs the adjusted limits pass.
          adjusted.input(Ext::vref) = published_vref;
          adjusted.input(Ext::pref) = published_pref;
          setState(adjusted.repca,
                   {{Vars::XQPI, 1.75},
                    {Vars::XQLAG, 0.0},
                    {Vars::XPPI, 2.25},
                    {Vars::PREF, 0.0}});
          success *= (adjusted.repca.evaluateResidual() == 0);
          success *= residualsMatch(adjusted.repca,
                                    {{Vars::XQLAG, 0.7}, {Vars::PREF, 4.5}},
                                    "adjusted command limits");
        }

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      /// Check every residual row against an independent numerical answer key.
      /// The expected values are literals, not a second implementation of REPCA.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        success *= (fixture.repca.evaluateResidual() == 0);

        // The state and inputs evaluate to V = 1, VLDC = 0.987, SFRZ = 1,
        // ERQ = 0.33, ERQLIM = 0.3, QPI = 0.05, EF = 0.2, EP = EPLIM = 0.05,
        // and PPI = 1.
        const std::array<VariableValue, Utilities::enum_size<Vars>()> expected_residuals{{
            {Vars::VMEAS, 1.235},
            {Vars::QMEAS, 0.45},
            {Vars::XQPI, 0.6},
            {Vars::XQLAG, 0.44},
            {Vars::PMEAS, 0.5},
            {Vars::XPPI, 0.69},
            {Vars::PREF, 0.6},
            {Vars::QEXT, -1.355},
            {Vars::PEXT, 0.05},
        }};

        success              *= (static_cast<size_t>(fixture.repca.getResidual().getSize())
                    == expected_residuals.size());
        const auto* residual  = fixture.repca.getResidual().getData();
        for (size_t row = 0; row < expected_residuals.size(); ++row)
        {
          const auto variable = expected_residuals[row].variable;
          if (index(variable) != row)
          {
            std::cout << "REPCA residual key position " << row << " names row "
                      << variableName(variable) << '\n';
            success = false;
          }
          success *= scalarMatches(residual[index(variable)],
                                   expected_residuals[row].value,
                                   variableName(variable));
        }

        return success.report(__func__);
      }

      /// Check reactive modes, smooth limits, antiwindup, and lead-lag behavior.
      TestOutcome reactiveControl()
      {
        TestStatus success = true;

        // The selected voltage reaches the filter row and the selected error
        // reaches the PI integrator row as Ki times the limited error.
        struct FlagCase
        {
          const char* label;
          bool        voltage_compensation;
          bool        voltage_reference;
          RealT       vmeas_rate;
          RealT       integrator_rate;
        };

        const std::array<FlagCase, 4> flag_cases{{
            {"droop/reactive-reference", false, false, 1.75, 1.11},
            {"droop/voltage-reference", false, true, 1.75, 0.51},
            {"line-drop/reactive-reference", true, false, 0.685, 1.11},
            {"line-drop/voltage-reference", true, true, 0.685, 0.51},
        }};
        for (const auto& test_case : flag_cases)
        {
          auto data                          = makeResidualData();
          data.parameters[Params::VcompFlag] = test_case.voltage_compensation;
          data.parameters[Params::RefFlag]   = test_case.voltage_reference;

          Fixture<ScalarT> fixture(data, kStateVr, kStateVi);
          fixture.attachAllInputs();
          setAnswerKeyInputs(fixture);
          fixture.input(Ext::vref) = 1.05;
          fixture.input(Ext::qref) = 0.25;

          success *= fixture.prepare(0.0, 0.0);
          setState(fixture.repca,
                   {{Vars::VMEAS, 0.85},
                    {Vars::QMEAS, 0.10},
                    {Vars::XQPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);

          const std::array<VariableValue, 2> expected_residuals{{
              {Vars::VMEAS, test_case.vmeas_rate},
              {Vars::XQPI, test_case.integrator_rate},
          }};
          success *= residualsMatch(fixture.repca,
                                    expected_residuals,
                                    test_case.label);
        }

        Fixture<ScalarT> fixture(makeResidualData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);

        // A voltage of zero must clear the freeze threshold by the same
        // margin the enabled probe clears it, so the threshold is raised. A
        // 0.2 limited error drives the integrator at 0.6 when enabled.
        {
          auto freeze_data                     = makeResidualData();
          freeze_data.parameters[Params::Vfrz] = 0.8;

          Fixture<ScalarT> freeze(freeze_data);
          freeze.attachAllInputs();
          setAnswerKeyInputs(freeze);
          success *= freeze.prepare(0.0, 0.0);

          const std::array<DrivenCase, 2> freeze_cases{{
              {0.0, 0.0},
              {1.6, 0.6},
          }};
          for (const auto& test_case : freeze_cases)
          {
            freeze.bus.Vr() = test_case.input;
            freeze.bus.Vi() = 0.0;
            freeze.bus.y().setDataUpdated();
            setState(freeze.repca, {{Vars::VMEAS, 0.82}, {Vars::XQPI, 0.0}});
            success *= (freeze.repca.evaluateResidual() == 0);
            success *= residualsMatch(freeze.repca,
                                      {{Vars::XQPI, test_case.expected}},
                                      "freeze gate");
          }
        }

        // The interior probe sits at the midpoint of the band, where the
        // smooth deadband cancels exactly. The deadbanded error passes the
        // error limit and reaches the integrator as Ki * erqdb.
        const std::array<DrivenCase, 3> deadband_cases{{
            {-0.27, -0.75},
            {0.005, 0.0},
            {0.28, 0.75},
        }};
        for (const auto& test_case : deadband_cases)
        {
          setState(fixture.repca,
                   {{Vars::VMEAS, 1.05 - test_case.input}, {Vars::XQPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XQPI, test_case.expected}},
                                    "reactive-power deadband");
        }

        // With Kp = 0 the PI output holds at XQPI, so the integrator reads Ki
        // times the limited error. Inputs are reactive errors beyond the band.
        {
          auto limit_data                   = makeResidualData();
          limit_data.parameters[Params::Kp] = 0.0;

          Fixture<ScalarT> limit(limit_data, kStateVr, kStateVi);
          limit.attachAllInputs();
          setAnswerKeyInputs(limit);
          success *= limit.prepare(0.0, 0.0);

          const std::array<DrivenCase, 3> error_limit_cases{{
              {-1.52, -2.1},
              {0.33, 0.9},
              {1.63, 2.4},
          }};
          for (const auto& test_case : error_limit_cases)
          {
            setState(limit.repca,
                     {{Vars::VMEAS, 1.05 - test_case.input}, {Vars::XQPI, 0.05}});
            success *= (limit.repca.evaluateResidual() == 0);
            success *= residualsMatch(limit.repca,
                                      {{Vars::XQPI, test_case.expected}},
                                      "reactive-power error limit");
          }
        }

        // With no error, the reactive lag reads the limited PI output.
        const std::array<DrivenCase, 3> command_limit_cases{{
            {-1.6, -0.32},
            {0.05, 0.02},
            {1.7, 0.36},
        }};
        for (const auto& test_case : command_limit_cases)
        {
          setState(fixture.repca,
                   {{Vars::VMEAS, 1.045},
                    {Vars::XQPI, test_case.input},
                    {Vars::XQLAG, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XQLAG, test_case.expected}},
                                    "reactive-power command limit");
        }

        // The integrator gate sees the unlimited PI input. Inside the limits
        // it passes the full rate; beyond a limit it passes restoring motion
        // and blocks outward motion. Cases give the unlimited PI input and the
        // limited error.
        const std::array<AntiWindupCase, 6> antiwindup_cases{{
            {-1.6, -0.4, 0.0},
            {-1.6, 0.4, 1.2},
            {0.05, -0.4, -1.2},
            {0.05, 0.4, 1.2},
            {1.7, -0.4, -1.2},
            {1.7, 0.4, 0.0},
        }};
        for (const auto& test_case : antiwindup_cases)
        {
          RealT band_edge = -0.02;
          if (test_case.error > 0.0)
          {
            band_edge = 0.03;
          }
          setState(fixture.repca,
                   {{Vars::VMEAS, 1.05 - (test_case.error + band_edge)},
                    {Vars::XQPI, test_case.output - 2.0 * test_case.error}});
          setDerivative(fixture.repca, {{Vars::XQPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XQPI, test_case.expected}},
                                    "reactive-power antiwindup");
        }

        setState(fixture.repca,
                 {{Vars::VMEAS, 1.045},
                  {Vars::XQPI, 0.27},
                  {Vars::XQLAG, 0.14},
                  {Vars::QEXT, 0.20}});
        setDerivative(fixture.repca, {{Vars::XQLAG, -0.04}});
        success *= (fixture.repca.evaluateResidual() == 0);
        const std::array<VariableValue, 2> lead_lag_residuals{{
            {Vars::XQLAG, 0.092},
            {Vars::QEXT, -0.624},
        }};
        success *= residualsMatch(fixture.repca,
                                  lead_lag_residuals,
                                  "reactive-command lead-lag");

        // A regulated voltage far below Vfrz freezes the integrator, leaving
        // it with no sensitivity to the error chain or the voltage.
        {
          Fixture<DependencyTracking::Variable> frozen(makeResidualData(), 0.1, 0.0);
          frozen.attachAllInputs();
          setAnswerKeyInputs(frozen);
          success *= frozen.prepare(0.0, 0.0);
          setState(frozen.repca, {{Vars::VMEAS, 0.72}, {Vars::XQPI, -0.55}});
          setDerivative(frozen.repca, {{Vars::XQPI, 0.0}});
          numberVariables(frozen);
          frozen.repca.updateTime(0.0, 1.0);
          success *= (frozen.repca.evaluateResidual() == 0);

          const DependencyTracking::Variable::DependencyMap expected{
              {2 * index(Vars::VMEAS), 0.0},               // @todo Remove these
              {2 * index(Vars::QMEAS), 0.0},               // @todo Remove these
              {2 * index(Vars::XQPI), 0.0},                // @todo Remove these
              {2 * index(Vars::XQPI) + 1, -1.0},           // @todo Remove these
              {2 * kBusVrColumn, 0.0},                     // @todo Remove these
              {2 * kBusViColumn, 0.0},                     // @todo Remove these
              {2 * externalColumn(index(Ext::vref)), 0.0}, // @todo Remove these
              {2 * externalColumn(index(Ext::qref)), 0.0}, // @todo Remove these
          };

          success *= jacobianRowMatches(
              frozen.repca.getResidual().getData()[index(Vars::XQPI)].getDependencies(),
              expected,
              index(Vars::XQPI),
              "frozen reactive-power integrator",
              kTol);
        }

        return success.report(__func__);
      }

      /// Check active-power modes, smooth limits, antiwindup, and command lag.
      TestOutcome activePowerControl()
      {
        TestStatus success = true;

        struct FlagCase
        {
          const char* label;
          bool        frequency_control;
          RealT       pext;
        };

        const std::array<FlagCase, 2> flag_cases{{
            {"disabled frequency control", false, -0.6},
            {"enabled frequency control", true, 0.2},
        }};
        for (const auto& test_case : flag_cases)
        {
          auto data                         = makeResidualData();
          data.parameters[Params::Freqflag] = test_case.frequency_control;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          setAnswerKeyInputs(fixture);
          success *= fixture.prepare(0.0, 0.0);
          setState(fixture.repca, {{Vars::PREF, 0.8}, {Vars::PEXT, 0.3}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::PEXT, test_case.pext}},
                                    test_case.label);
        }

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success                  *= fixture.prepare(0.0, 0.0);
        fixture.input(Ext::freq)  = 1.0;

        // The interior probe sits at the midpoint of the band, where the
        // smooth deadband cancels exactly. The deadbanded frequency error
        // passes the down (Ddn = 2) or up (Dup = 1) droop, adds to a 0.05
        // power error, and reaches the integrator as Kig times the result.
        const std::array<DrivenCase, 3> frequency_cases{{
            {-0.21, -0.63},
            {0.0025, 0.09},
            {0.315, 0.63},
        }};
        fixture.input(Ext::pref) = 0.175;
        for (const auto& test_case : frequency_cases)
        {
          fixture.input(Ext::freqref) = 1.0 + test_case.input;
          setState(fixture.repca, {{Vars::PMEAS, 0.3}, {Vars::XPPI, 1.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XPPI, test_case.expected}},
                                    "frequency deadband and droop");
        }

        // The remaining probes hold the frequency error at the deadband
        // midpoint and set the power error through the plant reference.
        fixture.input(Ext::freqref) = 1.0025;

        // With Kpg = 0 the PI output holds at XPPI, so the integrator reads
        // Kig times the limited power error.
        {
          auto limit_data                    = makeResidualData();
          limit_data.parameters[Params::Kpg] = 0.0;

          Fixture<ScalarT> limit(limit_data);
          limit.attachAllInputs();
          setAnswerKeyInputs(limit);
          success                   *= limit.prepare(0.0, 0.0);
          limit.input(Ext::freq)     = 1.0;
          limit.input(Ext::freqref)  = 1.0025;

          const std::array<DrivenCase, 3> error_limit_cases{{
              {-1.3, -0.9},
              {0.05, 0.09},
              {1.4, 1.08},
          }};
          for (const auto& test_case : error_limit_cases)
          {
            limit.input(Ext::pref) = 0.5 * (0.4 + test_case.input);
            setState(limit.repca, {{Vars::PMEAS, 0.4}, {Vars::XPPI, 1.0}});
            success *= (limit.repca.evaluateResidual() == 0);
            success *= residualsMatch(limit.repca,
                                      {{Vars::XPPI, test_case.expected}},
                                      "active-power error limit");
          }
        }

        // With no power error, the command lag reads the limited PI output.
        const std::array<DrivenCase, 3> command_limit_cases{{
            {-0.8, 0.0},
            {1.0, 2.0},
            {2.8, 4.0},
        }};
        fixture.input(Ext::pref) = 0.2;
        for (const auto& test_case : command_limit_cases)
        {
          setState(fixture.repca,
                   {{Vars::PMEAS, 0.4},
                    {Vars::XPPI, test_case.input},
                    {Vars::PREF, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::PREF, test_case.expected}},
                                    "active-power command limit");
        }

        // The integrator gate sees the unlimited PI input: full rate inside
        // the limits, restoring motion passed and outward motion blocked
        // beyond them. Cases give the unlimited PI input and the limited power
        // error.
        const std::array<AntiWindupCase, 6> antiwindup_cases{{
            {-0.8, -0.3, 0.0},
            {-0.8, 0.3, 0.54},
            {1.0, -0.3, -0.54},
            {1.0, 0.3, 0.54},
            {2.8, -0.3, -0.54},
            {2.8, 0.3, 0.0},
        }};
        for (const auto& test_case : antiwindup_cases)
        {
          fixture.input(Ext::pref) = 0.5 * (0.4 + test_case.error);
          setState(fixture.repca,
                   {{Vars::PMEAS, 0.4},
                    {Vars::XPPI, test_case.output - 1.7 * test_case.error}});
          setDerivative(fixture.repca, {{Vars::XPPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XPPI, test_case.expected}},
                                    "active-power antiwindup");
        }

        fixture.input(Ext::pref) = 0.2;
        setState(fixture.repca,
                 {{Vars::PMEAS, 0.4}, {Vars::XPPI, 0.66}, {Vars::PREF, 0.60}});
        setDerivative(fixture.repca, {{Vars::PREF, 0.05}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{Vars::PREF, 0.07}},
                                  "active-power command lag");

        return success.report(__func__);
      }

      /// Check every differential residual with nonzero explicit derivatives.
      TestOutcome derivatives()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        setDerivative(fixture.repca,
                      {{Vars::VMEAS, 0.2},
                       {Vars::QMEAS, -0.1},
                       {Vars::XQPI, 0.4},
                       {Vars::XQLAG, -0.3},
                       {Vars::PMEAS, 0.6},
                       {Vars::XPPI, -0.5},
                       {Vars::PREF, 0.8}});
        success *= (fixture.repca.evaluateResidual() == 0);
        const std::array<VariableValue, 7> expected_residuals{{
            {Vars::VMEAS, 1.135},
            {Vars::QMEAS, 0.35},
            {Vars::XQPI, 0.5},
            {Vars::XQLAG, 0.34},
            {Vars::PMEAS, 0.4},
            {Vars::XPPI, 0.59},
            {Vars::PREF, 0.5},
        }};
        success *= residualsMatch(fixture.repca,
                                  expected_residuals,
                                  "explicit derivatives");

        return success.report(__func__);
      }

      /// Check every dependency-tracking Jacobian row against an independent
      /// numerical and structural answer key, with both selector settings and
      /// a non-unit alpha.
      TestOutcome dependencyTracking()
      {
        TestStatus success = true;

        const auto data       = makeResidualData();
        const auto dependency = dependencyTrackingJacobian(data, success);

        success *= jacobianMatches(dependency,
                                   expectedJacobian(),
                                   "dependency tracking",
                                   kTol);

        auto all_flags_off_data                          = data;
        all_flags_off_data.parameters[Params::VcompFlag] = false;
        all_flags_off_data.parameters[Params::RefFlag]   = false;
        all_flags_off_data.parameters[Params::Freqflag]  = false;
        const auto all_flags_off =
            dependencyTrackingJacobian(all_flags_off_data, success);
        success *= jacobianMatches(all_flags_off,
                                   expectedJacobianAllFlagsOff(),
                                   "all-flags-off dependency tracking",
                                   kTol);

        const auto nonunit_alpha_dependency =
            dependencyTrackingJacobian(data, success, kNonunitAlpha);
        success *= jacobianMatches(nonunit_alpha_dependency,
                                   expectedJacobianNonunitAlpha(),
                                   "non-unit-alpha dependency tracking",
                                   kTol);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// One rich state, both selector settings, and a non-unit alpha drive both
      /// sensitivity paths; every Enzyme CSR row must match dependency tracking.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto data = makeResidualData();

        success *= jacobianMatches(enzymeJacobian(data, success),
                                   dependencyTrackingJacobian(data, success),
                                   "Enzyme versus dependency tracking",
                                   kTol);

        auto all_flags_off_data                           = data;
        all_flags_off_data.parameters[Params::VcompFlag]  = false;
        all_flags_off_data.parameters[Params::RefFlag]    = false;
        all_flags_off_data.parameters[Params::Freqflag]   = false;
        success                                          *= jacobianMatches(
            enzymeJacobian(all_flags_off_data, success),
            dependencyTrackingJacobian(all_flags_off_data, success),
            "all-flags-off Enzyme versus dependency tracking",
            kTol);

        success *= jacobianMatches(
            enzymeJacobian(data, success, kNonunitAlpha),
            dependencyTrackingJacobian(data, success, kNonunitAlpha),
            "non-unit-alpha Enzyme versus dependency tracking",
            kTol);

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Controller::RepcaParameters;
      using Vars   = PhasorDynamics::Controller::RepcaInternalVariables;
      using Mon    = PhasorDynamics::Controller::RepcaMonitorableVariables;
      using Data   = PhasorDynamics::Controller::RepcaData<RealT, IdxT>;
      using Ext    = typename Data::SignalInputs;
      using RepcaT = PhasorDynamics::Controller::Repca<ScalarT, IdxT>;

      static constexpr size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(Ext variable)
      {
        return static_cast<size_t>(variable);
      }

      struct VariableValue
      {
        Vars  variable;
        RealT value;
      };

      struct DrivenCase
      {
        RealT input;
        RealT expected;
      };

      struct AntiWindupCase
      {
        RealT output;
        RealT error;
        RealT expected;
      };

      enum class NonfiniteTarget
      {
        NONE,
        INPUT,
        BUS_VOLTAGE
      };

      /// Owns the regulated bus, REPCA, assigned command nodes, and attached
      /// input nodes. Signal storage precedes the model so every referenced node
      /// outlives REPCA; copying would invalidate the model and node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, Utilities::enum_size<Ext>()>                                   input_values_{};
        std::array<IdxT, Utilities::enum_size<Ext>()>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, Utilities::enum_size<Ext>()> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> qext_node_;
        PhasorDynamics::SignalNode<T, IdxT> pext_node_;

      public:
        explicit Fixture(const Data& data,
                         RealT       vr                     = 1.0,
                         RealT       vi                     = 0.0,
                         RealT       system_va_base         = 100.0e6,
                         bool        assign_command_outputs = true)
          : bus(static_cast<T>(vr), static_cast<T>(vi)),
            repca(&bus, data)
        {
          repca.setSystemBase(60.0, system_va_base);
          if (assign_command_outputs)
          {
            repca.getPorts().out.template port<Data::SignalOutputs::qext>().connect(&qext_node_);
            repca.getPorts().out.template port<Data::SignalOutputs::pext>().connect(&pext_node_);
          }
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachRequiredInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = repca.size() + bus.size();
          for (auto variant : Utilities::enum_values<Ext>())
          {
            const auto port      = static_cast<size_t>(variant);
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].link(&input_values_[port], &input_indices_[port]);
          }

          repca.getPorts().in.template port<Data::SignalInputs::ir>().connect(&input_nodes_[index(Ext::ir)]);
          repca.getPorts().in.template port<Data::SignalInputs::ii>().connect(&input_nodes_[index(Ext::ii)]);
          repca.getPorts().in.template port<Data::SignalInputs::p>().connect(&input_nodes_[index(Ext::p)]);
          repca.getPorts().in.template port<Data::SignalInputs::q>().connect(&input_nodes_[index(Ext::q)]);
        }

        void attachAllInputs(RealT initial_value    = 0.0,
                             bool  attach_frequency = true)
        {
          attachRequiredInputs(initial_value);

          if (attach_frequency)
          {
            repca.getPorts().in.template port<Data::SignalInputs::freq>().connect(&input_nodes_[index(Ext::freq)]);
          }
          repca.getPorts().in.template port<Data::SignalInputs::vref>().connect(&input_nodes_[index(Ext::vref)]);
          repca.getPorts().in.template port<Data::SignalInputs::pref>().connect(&input_nodes_[index(Ext::pref)]);
          repca.getPorts().in.template port<Data::SignalInputs::qref>().connect(&input_nodes_[index(Ext::qref)]);
          repca.getPorts().in.template port<Data::SignalInputs::freqref>().connect(&input_nodes_[index(Ext::freqref)]);
        }

        void setCommands(RealT qext, RealT pext)
        {
          auto* y              = repca.y().getData();
          y[index(Vars::QEXT)] = static_cast<T>(qext);
          y[index(Vars::PEXT)] = static_cast<T>(pext);
          repca.y().setDataUpdated();
        }

        /// Arrange the allocation, verification, bus, and command prerequisites.
        bool prepare(RealT qext, RealT pext)
        {
          const bool success = (bus.allocate() == 0) && (repca.allocate() == 0)
                               && (repca.verify() == 0) && (bus.initialize() == 0);
          if (!success)
          {
            std::cout << "REPCA fixture preparation failed\n";
            return false;
          }
          setCommands(qext, pext);
          return true;
        }

        bool initialize(RealT qext, RealT pext)
        {
          if (!prepare(qext, pext))
          {
            return false;
          }
          if (repca.initialize() != 0)
          {
            std::cout << "REPCA initialization failed\n";
            return false;
          }
          return true;
        }

        T qext() const
        {
          return repca.y().getData()[index(Vars::QEXT)];
        }

        T pext() const
        {
          return repca.y().getData()[index(Vars::PEXT)];
        }

        T& input(Ext port)
        {
          return input_values_[index(port)];
        }

        IdxT inputIndex(Ext port) const
        {
          return input_indices_[index(port)];
        }

        PhasorDynamics::Bus<T, IdxT>               bus;
        PhasorDynamics::Controller::Repca<T, IdxT> repca;
      };

      static constexpr RealT kStateVr      = 0.6;
      static constexpr RealT kStateVi      = 0.8;
      static constexpr RealT kNonunitAlpha = 2.5;

      static constexpr size_t kBusVrColumn        = Utilities::enum_size<Vars>();
      static constexpr size_t kBusViColumn        = kBusVrColumn + 1;
      static constexpr size_t kExternalColumnBase = kBusViColumn + 1;

      static constexpr size_t externalColumn(size_t port)
      {
        return kExternalColumnBase + port;
      }

      Data makeMinimalData() const
      {
        Data data;
        data.device_class          = "Repca";
        data.disambiguation_string = "repca_test";
        data.monitored_variables.insert(Mon::qext);
        data.monitored_variables.insert(Mon::pext);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::qmeas);
        data.monitored_variables.insert(Mon::pmeas);
        return data;
      }

      Data makeExplicitDefaultData() const
      {
        auto data                          = makeMinimalData();
        data.parameters[Params::mva]       = 100.0;
        data.parameters[Params::VcompFlag] = true;
        data.parameters[Params::RefFlag]   = true;
        data.parameters[Params::Freqflag]  = false;
        data.parameters[Params::Tfltr]     = 0.05;
        data.parameters[Params::Vfrz]      = 0.7;
        data.parameters[Params::Rc]        = 0.0;
        data.parameters[Params::Xc]        = 0.0;
        data.parameters[Params::Kc]        = 1.0;
        data.parameters[Params::dbdlow]    = 0.0;
        data.parameters[Params::dbdupper]  = 0.0;
        data.parameters[Params::emax]      = 1.0;
        data.parameters[Params::emin]      = -1.0;
        data.parameters[Params::Kp]        = 10.0;
        data.parameters[Params::Ki]        = 10.0;
        data.parameters[Params::Qmax]      = 1.0;
        data.parameters[Params::Qmin]      = -1.0;
        data.parameters[Params::Tft]       = 0.0;
        data.parameters[Params::Tfv]       = 3.0;
        data.parameters[Params::Tp]        = 0.0;
        data.parameters[Params::fdbd1]     = 0.0;
        data.parameters[Params::fdbd2]     = 0.0;
        data.parameters[Params::Ddn]       = 20.0;
        data.parameters[Params::Dup]       = 0.0;
        data.parameters[Params::femax]     = 1.0;
        data.parameters[Params::femin]     = -1.0;
        data.parameters[Params::Kpg]       = 10.0;
        data.parameters[Params::Kig]       = 10.0;
        data.parameters[Params::Pmax]      = 2.0;
        data.parameters[Params::Pmin]      = 0.0;
        data.parameters[Params::Tlag]      = 3.0;
        return data;
      }

      Data makeData() const
      {
        auto data                         = makeExplicitDefaultData();
        data.parameters[Params::Freqflag] = true;
        data.parameters[Params::Tp]       = 0.05;
        return data;
      }

      /// Distinct nonzero values for every parameter. The limiter bands are
      /// wide enough, and the lag reciprocals exact enough, for probe states
      /// to clear every smooth transition on an exact decimal.
      Data makeResidualData() const
      {
        auto data                         = makeData();
        data.parameters[Params::mva]      = 50.0;
        data.parameters[Params::Tfltr]    = 0.2;
        data.parameters[Params::Rc]       = 0.02;
        data.parameters[Params::Xc]       = 0.03;
        data.parameters[Params::Kc]       = 0.4;
        data.parameters[Params::dbdlow]   = -0.02;
        data.parameters[Params::dbdupper] = 0.03;
        data.parameters[Params::emax]     = 0.8;
        data.parameters[Params::emin]     = -0.7;
        data.parameters[Params::Kp]       = 2.0;
        data.parameters[Params::Ki]       = 3.0;
        data.parameters[Params::Qmax]     = 0.9;
        data.parameters[Params::Qmin]     = -0.8;
        data.parameters[Params::Tft]      = 0.2;
        data.parameters[Params::Tfv]      = 2.5;
        data.parameters[Params::Tp]       = 0.4;
        data.parameters[Params::fdbd1]    = -0.01;
        data.parameters[Params::fdbd2]    = 0.015;
        data.parameters[Params::Ddn]      = 2.0;
        data.parameters[Params::Dup]      = 1.0;
        data.parameters[Params::femax]    = 0.6;
        data.parameters[Params::femin]    = -0.5;
        data.parameters[Params::Kpg]      = 1.7;
        data.parameters[Params::Kig]      = 1.8;
        data.parameters[Params::Pmax]     = 2.0;
        data.parameters[Params::Tlag]     = 0.5;
        return data;
      }

      /// Both deadbands and both error limits are symmetric and the droop
      /// gains are equal, so an operating point with no error reconstructs
      /// exactly; the commands and the freeze threshold clear their limits.
      Data makeInitializationData() const
      {
        auto data                         = makeResidualData();
        data.parameters[Params::Vfrz]     = 0.2;
        data.parameters[Params::dbdupper] = 0.02;
        data.parameters[Params::emin]     = -0.8;
        data.parameters[Params::Qmax]     = 1.5;
        data.parameters[Params::fdbd1]    = -0.015;
        data.parameters[Params::Dup]      = 2.0;
        data.parameters[Params::femin]    = -0.6;
        return data;
      }

      template <typename T>
      void setInitializationInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::ir)   = static_cast<T>(0.2);
        fixture.input(Ext::ii)   = static_cast<T>(-0.1);
        fixture.input(Ext::p)    = static_cast<T>(0.4);
        fixture.input(Ext::q)    = static_cast<T>(0.1);
        fixture.input(Ext::freq) = static_cast<T>(0.99);
      }

      /// On the unit terminal (0.6, 0.8), the branch current drops the
      /// compensated voltage along the terminal phasor to 0.987, and both the
      /// voltage and reactive references sit 0.33 above their measurements.
      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::ir)      = static_cast<T>(0.18);
        fixture.input(Ext::ii)      = static_cast<T>(-0.01);
        fixture.input(Ext::p)       = static_cast<T>(0.35);
        fixture.input(Ext::q)       = static_cast<T>(0.25);
        fixture.input(Ext::freq)    = static_cast<T>(1.0);
        fixture.input(Ext::vref)    = static_cast<T>(1.05);
        fixture.input(Ext::pref)    = static_cast<T>(0.075);
        fixture.input(Ext::qref)    = static_cast<T>(0.39);
        fixture.input(Ext::freqref) = static_cast<T>(1.215);
      }

      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Controller::Repca<T, IdxT>& repca) const
      {
        // Every smooth-transition argument keeps a saturation margin, and
        // every clamp that must pass its input through sits at the midpoint
        // of its limits, so each row carries its ideal value.
        setState(repca,
                 {{Vars::VMEAS, 0.72},
                  {Vars::QMEAS, 0.45},
                  {Vars::XQPI, -0.55},
                  {Vars::XQLAG, -0.05},
                  {Vars::PMEAS, 0.3},
                  {Vars::XPPI, 0.915},
                  {Vars::PREF, 0.35},
                  {Vars::QEXT, 0.25},
                  {Vars::PEXT, 0.15}});
        setDerivative(repca,
                      {{Vars::VMEAS, 0.1},
                       {Vars::QMEAS, -0.2},
                       {Vars::XQPI, 0.3},
                       {Vars::XQLAG, -0.4},
                       {Vars::PMEAS, 0.5},
                       {Vars::XPPI, -0.6},
                       {Vars::PREF, 0.7}});
      }

      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), 0.9, 0.4);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), 0.9, 0.4);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        implicit_defaults.input(Ext::p)    = 0.2;
        implicit_defaults.input(Ext::q)    = 0.1;
        implicit_defaults.input(Ext::freq) = 1.0;
        explicit_defaults.input(Ext::p)    = 0.2;
        explicit_defaults.input(Ext::q)    = 0.1;
        explicit_defaults.input(Ext::freq) = 1.0;

        bool success = implicit_defaults.initialize(0.1, 0.2)
                       && explicit_defaults.initialize(0.1, 0.2);
        if (!success)
        {
          std::cout << "REPCA documented-default comparison failed to initialize\n";
          return false;
        }

        if (implicit_defaults.repca.evaluateResidual() != 0)
        {
          success = false;
        }
        if (explicit_defaults.repca.evaluateResidual() != 0)
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.repca.y(),
                          explicit_defaults.repca.y(),
                          "documented-default state"))
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.repca.yp(),
                          explicit_defaults.repca.yp(),
                          "documented-default derivative"))
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.repca.getResidual(),
                          explicit_defaults.repca.getResidual(),
                          "documented-default residual"))
        {
          success = false;
        }
        for (auto variant : Utilities::enum_values<Ext>())
        {
          const auto port     = static_cast<size_t>(variant);
          const auto variable = static_cast<Ext>(port);
          if (!rowMatches(implicit_defaults.input(variable),
                          explicit_defaults.input(variable),
                          "documented-default signal",
                          port,
                          ""))
          {
            success = false;
          }
        }

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.repca);
        setAnswerKeyState(explicit_defaults.repca);
        if (implicit_defaults.repca.evaluateResidual() != 0)
        {
          success = false;
        }
        if (explicit_defaults.repca.evaluateResidual() != 0)
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.repca.getResidual(),
                          explicit_defaults.repca.getResidual(),
                          "documented-default dynamic residual"))
        {
          success = false;
        }
        return success;
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        return fixture.repca.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected() const
      {
        Fixture<ScalarT> fixture(makeData());
        fixture.attachAllInputs();
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        fixture.repca.getPorts().in.template port<variable>().connect(&unlinked_node);
        return fixture.repca.verify() > 0;
      }

      /// Fill state and derivative with a recognizable ramp, restoring the
      /// aliased commands, so any write by a rejected initialization shows.
      void poisonState(Fixture<ScalarT>& fixture, RealT qext, RealT pext) const
      {
        auto* y  = fixture.repca.y().getData();
        auto* yp = fixture.repca.yp().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          y[row]  = 0.125 + 0.01 * static_cast<RealT>(row);
          yp[row] = -0.25 - 0.01 * static_cast<RealT>(row);
        }
        fixture.setCommands(qext, pext);
        fixture.repca.yp().setDataUpdated();
      }

      bool initializationRejectedAtomically(const Data&     data,
                                            RealT           qext,
                                            RealT           pext,
                                            const char*     label,
                                            NonfiniteTarget target        = NonfiniteTarget::NONE,
                                            RealT           initial_vr    = 0.8,
                                            RealT           initial_vi    = 0.6,
                                            Ext             poisoned_port = Ext::freq,
                                            RealT           poison_value =
                                                std::numeric_limits<RealT>::infinity()) const
      {
        Fixture<ScalarT> fixture(data, initial_vr, initial_vi);
        fixture.attachAllInputs(77.0);
        setInitializationInputs(fixture);
        if (!fixture.prepare(qext, pext))
        {
          return false;
        }

        if (target == NonfiniteTarget::INPUT)
        {
          fixture.input(poisoned_port) = poison_value;
        }
        if (target == NonfiniteTarget::BUS_VOLTAGE)
        {
          fixture.bus.Vr() = poison_value;
          fixture.bus.y().setDataUpdated();
        }

        poisonState(fixture, qext, pext);

        const auto                                     y_before   = copyVector(fixture.repca.y());
        const auto                                     yp_before  = copyVector(fixture.repca.yp());
        const auto                                     bus_before = copyVector(fixture.bus.y());
        std::array<RealT, Utilities::enum_size<Ext>()> inputs_before{};
        for (auto variant : Utilities::enum_values<Ext>())
        {
          const auto port     = static_cast<size_t>(variant);
          inputs_before[port] = fixture.input(static_cast<Ext>(port));
        }

        bool success = true;
        if (fixture.repca.initialize() == 0)
        {
          std::cout << "Expected REPCA initialization rejection: " << label << '\n';
          success = false;
        }

        if (!scalarPreserved(fixture.qext(), qext, "rejected qext preservation"))
        {
          success = false;
        }
        if (!scalarPreserved(fixture.pext(), pext, "rejected pext preservation"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.repca.y(), y_before, "state"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.repca.yp(), yp_before, "derivative"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.bus.y(), bus_before, "bus state"))
        {
          success = false;
        }
        for (auto variant : Utilities::enum_values<Ext>())
        {
          const auto port = static_cast<size_t>(variant);
          if (!valueUnchanged(fixture.input(static_cast<Ext>(port)),
                              inputs_before[port],
                              "external signal",
                              port))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename T>
      void setState(PhasorDynamics::Controller::Repca<T, IdxT>& repca,
                    std::initializer_list<VariableValue>        values) const
      {
        auto* y = repca.y().getData();
        for (const auto& [variable, value] : values)
        {
          y[index(variable)] = static_cast<T>(value);
        }
        repca.y().setDataUpdated();
      }

      template <typename T>
      void setDerivative(PhasorDynamics::Controller::Repca<T, IdxT>& repca,
                         std::initializer_list<VariableValue>        values) const
      {
        auto* yp = repca.yp().getData();
        for (const auto& [variable, value] : values)
        {
          yp[index(variable)] = static_cast<T>(value);
        }
        repca.yp().setDataUpdated();
      }

      static const char* variableName(Vars variable)
      {
        static constexpr std::array<const char*, Utilities::enum_size<Vars>()> names{{
            "VMEAS",
            "QMEAS",
            "XQPI",
            "XQLAG",
            "PMEAS",
            "XPPI",
            "PREF",
            "QEXT",
            "PEXT",
        }};
        return names[index(variable)];
      }

      static bool variableMatches(RealT       actual,
                                  RealT       expected,
                                  const char* what,
                                  Vars        variable,
                                  const char* context,
                                  RealT       tolerance = kTol)
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "REPCA " << what << ' ' << variableName(variable);
        if (context[0] != '\0')
        {
          std::cout << ' ' << context;
        }
        std::cout << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      static bool rowMatches(RealT       actual,
                             RealT       expected,
                             const char* what,
                             size_t      row,
                             const char* context,
                             RealT       tolerance = kTol)
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "REPCA " << what << " row " << row;
        if (context[0] != '\0')
        {
          std::cout << ' ' << context;
        }
        std::cout << " mismatch: " << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual
                  << " != " << expected << '\n';
        return false;
      }

      bool scalarMatches(RealT       actual,
                         RealT       expected,
                         const char* label,
                         RealT       tolerance = kTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << label << " mismatch: " << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual
                  << " != " << expected << '\n';
        return false;
      }

      bool monitorMatches(const RepcaT&               repca,
                          const std::array<RealT, 5>& expected,
                          const char*                 context) const
      {
        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(repca.getMonitor());
        std::stringstream output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string header;
        std::string values_line;
        std::getline(output, header);
        std::getline(output, values_line);

        bool success =
            header == "t,Repca_repca_test_qext,Repca_repca_test_pext,"
                      "Repca_repca_test_vmeas,Repca_repca_test_qmeas,"
                      "Repca_repca_test_pmeas";

        const auto values = Tokenizer<RealT>(values_line, ',')();
        if (values.size() != expected.size() + 1)
        {
          std::cout << "REPCA monitor emitted " << values.size()
                    << " values instead of " << expected.size() + 1 << '\n';
          return false;
        }

        for (size_t i = 0; i < expected.size(); ++i)
        {
          if (!rowMatches(values[i + 1],
                          expected[i],
                          "monitor",
                          i,
                          context))
          {
            success = false;
          }
        }
        return success;
      }

      /// A value retains exactly what its owner supplied, including signed
      /// infinities and NaN.
      static bool preserved(RealT actual, RealT expected)
      {
        if (std::isnan(expected))
        {
          return std::isnan(actual);
        }
        return actual == expected;
      }

      bool scalarPreserved(RealT actual, RealT expected, const char* label) const
      {
        if (preserved(actual, expected))
        {
          return true;
        }
        std::cout << label << " changed: " << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual
                  << " != " << expected << '\n';
        return false;
      }

      static bool valueUnchanged(RealT       actual,
                                 RealT       expected,
                                 const char* what,
                                 size_t      index)
      {
        if (preserved(actual, expected))
        {
          return true;
        }
        std::cout << "REPCA " << what << ' ' << index
                  << " changed: " << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual
                  << " != " << expected << '\n';
        return false;
      }

      template <typename VectorT, typename ValuesT>
      bool rowsMatch(const VectorT& vector,
                     const ValuesT& values,
                     const char*    what,
                     const char*    context) const
      {
        bool        success       = true;
        const auto* vector_values = vector.getData();
        for (const auto& [variable, expected] : values)
        {
          const size_t row = index(variable);

          if (!variableMatches(static_cast<RealT>(vector_values[row]),
                               expected,
                               what,
                               variable,
                               context))
          {
            success = false;
          }
        }
        return success;
      }

      bool residualsMatch(const RepcaT&                        repca,
                          std::initializer_list<VariableValue> values,
                          const char*                          context = "") const
      {
        return rowsMatch(repca.getResidual(), values, "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const RepcaT&                          repca,
                          const std::array<VariableValue, size>& values,
                          const char*                            context = "") const
      {
        return rowsMatch(repca.getResidual(), values, "residual", context);
      }

      bool stateMatches(const RepcaT&                        repca,
                        std::initializer_list<VariableValue> values,
                        const char*                          context = "") const
      {
        return rowsMatch(repca.y(), values, "state", context);
      }

      template <size_t size>
      bool stateMatches(const RepcaT&                          repca,
                        const std::array<VariableValue, size>& values,
                        const char*                            context = "") const
      {
        return rowsMatch(repca.y(), values, "state", context);
      }

      bool allResidualsWithinInitTolerance(const RepcaT& repca) const
      {
        bool        success = true;
        const auto* f       = repca.getResidual().getData();
        const auto* yp      = repca.yp().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          const auto variable = static_cast<Vars>(row);
          if (!variableMatches(f[row],
                               0.0,
                               "residual",
                               variable,
                               "at rest",
                               RepcaT::INITIALIZATION_TOLERANCE))
          {
            success = false;
          }
          if (!valueUnchanged(yp[row], 0.0, "derivative", row))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename VectorT>
      std::vector<RealT> copyVector(const VectorT& vector) const
      {
        const auto*        values = vector.getData();
        std::vector<RealT> snapshot(static_cast<size_t>(vector.getSize()));
        for (size_t row = 0; row < snapshot.size(); ++row)
        {
          snapshot[row] = static_cast<RealT>(values[row]);
        }
        return snapshot;
      }

      template <typename VectorT>
      bool vectorUnchanged(const VectorT&            vector,
                           const std::vector<RealT>& snapshot,
                           const char*               what) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t row = 0; row < snapshot.size(); ++row)
        {
          if (!valueUnchanged(static_cast<RealT>(values[row]),
                              snapshot[row],
                              what,
                              row))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename LeftVectorT, typename RightVectorT>
      bool vectorsMatch(const LeftVectorT&  left,
                        const RightVectorT& right,
                        const char*         what) const
      {
        if (left.getSize() != right.getSize())
        {
          std::cout << "REPCA " << what << " size mismatch\n";
          return false;
        }
        bool        success      = true;
        const auto* left_values  = left.getData();
        const auto* right_values = right.getData();
        for (size_t row = 0; row < static_cast<size_t>(left.getSize()); ++row)
        {
          if (!rowMatches(static_cast<RealT>(left_values[row]),
                          static_cast<RealT>(right_values[row]),
                          what,
                          row,
                          ""))
          {
            success = false;
          }
        }
        return success;
      }

      /// Each row includes the structural entries of every input its
      /// evaluated chain reads, with zero value where a mode mask or a
      /// saturated smooth gate removes the sensitivity.
      std::vector<DependencyTracking::Variable::DependencyMap> expectedJacobian() const
      {
        return {
            {{index(Vars::VMEAS), -6.0},
             {kBusVrColumn, 3.0},
             {kBusViColumn, 4.0},
             {externalColumn(index(Ext::ir)), -0.36},
             {externalColumn(index(Ext::ii)), 0.02},
             {externalColumn(index(Ext::q)), 0.0}},
            {{index(Vars::QMEAS), -6.0}, {externalColumn(index(Ext::q)), 10.0}},
            {{index(Vars::VMEAS), -3.0},
             {index(Vars::QMEAS), 0.0},
             {index(Vars::XQPI), -1.0},
             {kBusVrColumn, 0.0},
             {kBusViColumn, 0.0},
             {externalColumn(index(Ext::vref)), 3.0},
             {externalColumn(index(Ext::qref)), 0.0}},
            {{index(Vars::VMEAS), -0.8},
             {index(Vars::QMEAS), 0.0},
             {index(Vars::XQPI), 0.4},
             {index(Vars::XQLAG), -1.4},
             {externalColumn(index(Ext::vref)), 0.8},
             {externalColumn(index(Ext::qref)), 0.0}},
            {{index(Vars::PMEAS), -3.5}, {externalColumn(index(Ext::p)), 5.0}},
            {{index(Vars::PMEAS), -1.8},
             {index(Vars::XPPI), -1.0},
             {externalColumn(index(Ext::freq)), -1.8},
             {externalColumn(index(Ext::pref)), 3.6},
             {externalColumn(index(Ext::freqref)), 1.8}},
            {{index(Vars::PMEAS), -3.4},
             {index(Vars::XPPI), 2.0},
             {index(Vars::PREF), -3.0},
             {externalColumn(index(Ext::freq)), -3.4},
             {externalColumn(index(Ext::pref)), 6.8},
             {externalColumn(index(Ext::freqref)), 3.4}},
            {{index(Vars::VMEAS), -0.4},
             {index(Vars::QMEAS), 0.0},
             {index(Vars::XQPI), 0.2},
             {index(Vars::XQLAG), 2.3},
             {index(Vars::QEXT), -5.0},
             {externalColumn(index(Ext::vref)), 0.4},
             {externalColumn(index(Ext::qref)), 0.0}},
            {{index(Vars::PREF), 1.0}, {index(Vars::PEXT), -2.0}},
        };
      }

      std::vector<DependencyTracking::Variable::DependencyMap> expectedJacobianAllFlagsOff() const
      {
        auto expected                = expectedJacobian();
        expected[index(Vars::VMEAS)] = {
            {index(Vars::VMEAS), -6.0},
            {kBusVrColumn, 3.0},
            {kBusViColumn, 4.0},
            {externalColumn(index(Ext::ir)), 0.0},
            {externalColumn(index(Ext::ii)), 0.0},
            {externalColumn(index(Ext::q)), 4.0},
        };
        expected[index(Vars::XQPI)] = {
            {index(Vars::VMEAS), 0.0},
            {index(Vars::QMEAS), -3.0},
            {index(Vars::XQPI), -1.0},
            {kBusVrColumn, 0.0},
            {kBusViColumn, 0.0},
            {externalColumn(index(Ext::vref)), 0.0},
            {externalColumn(index(Ext::qref)), 6.0},
        };
        expected[index(Vars::XQLAG)] = {
            {index(Vars::VMEAS), 0.0},
            {index(Vars::QMEAS), -0.8},
            {index(Vars::XQPI), 0.4},
            {index(Vars::XQLAG), -1.4},
            {externalColumn(index(Ext::vref)), 0.0},
            {externalColumn(index(Ext::qref)), 1.6},
        };
        expected[index(Vars::QEXT)] = {
            {index(Vars::VMEAS), 0.0},
            {index(Vars::QMEAS), -0.4},
            {index(Vars::XQPI), 0.2},
            {index(Vars::XQLAG), 2.3},
            {index(Vars::QEXT), -5.0},
            {externalColumn(index(Ext::vref)), 0.0},
            {externalColumn(index(Ext::qref)), 0.8},
        };
        expected[index(Vars::PEXT)] = {
            {index(Vars::PREF), 0.0},
            {index(Vars::PEXT), -2.0},
        };
        return expected;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> expectedJacobianNonunitAlpha() const
      {
        auto expected                                    = expectedJacobian();
        expected[index(Vars::VMEAS)][index(Vars::VMEAS)] = -7.5;
        expected[index(Vars::QMEAS)][index(Vars::QMEAS)] = -7.5;
        expected[index(Vars::XQPI)][index(Vars::XQPI)]   = -2.5;
        expected[index(Vars::XQLAG)][index(Vars::XQLAG)] = -2.9;
        expected[index(Vars::PMEAS)][index(Vars::PMEAS)] = -5.0;
        expected[index(Vars::XPPI)][index(Vars::XPPI)]   = -2.5;
        expected[index(Vars::PREF)][index(Vars::PREF)]   = -4.5;
        return expected;
      }

      bool jacobianRowMatches(
          const DependencyTracking::Variable::DependencyMap& actual,
          const DependencyTracking::Variable::DependencyMap& expected,
          size_t                                             row,
          const char*                                        source,
          RealT                                              tolerance) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }

        std::cout << "REPCA " << source << " Jacobian row " << row
                  << " mismatch\n";
        return false;
      }

      bool jacobianMatches(
          const std::vector<DependencyTracking::Variable::DependencyMap>& actual,
          const std::vector<DependencyTracking::Variable::DependencyMap>& expected,
          const char*                                                     source,
          RealT                                                           tolerance) const
      {
        if (actual.size() != expected.size())
        {
          std::cout << "REPCA " << source << " Jacobian row-count mismatch\n";
          return false;
        }

        bool success = true;
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          if (!jacobianRowMatches(actual[row], expected[row], row, source, tolerance))
          {
            success = false;
          }
        }
        return success;
      }

      /// The reactive error the published reference reproduces at rest, on
      /// the component base (makeInitializationData() halves the power base).
      RealT reactiveError(Fixture<ScalarT>& fixture, bool voltage_reference) const
      {
        const auto* y = fixture.repca.y().getData();
        if (voltage_reference)
        {
          return fixture.input(Ext::vref) - y[index(Vars::VMEAS)];
        }
        return 2.0 * fixture.input(Ext::qref) - y[index(Vars::QMEAS)];
      }

      /// @todo Remove and setup the test to not rely on explicit variable numbering
      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y     = fixture.repca.y().getData();
        auto* yp    = fixture.repca.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          y[row].setVariableNumber(2 * row);
          yp[row].setVariableNumber(2 * row + 1);
        }
        for (size_t row = 0; row < static_cast<size_t>(fixture.bus.size()); ++row)
        {
          bus_y[row].setVariableNumber(2 * (kBusVrColumn + row));
        }
        for (auto variable : Utilities::enum_values<Ext>())
        {
          fixture.input(variable).setVariableNumber(2 * fixture.inputIndex(variable));
        }

        fixture.repca.y().setDataUpdated();
        fixture.repca.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          const Data& data,
          TestStatus& success,
          RealT       alpha = 1.0) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        numberVariables(fixture);
        fixture.repca.updateTime(0.0, alpha);
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= (fixture.repca.evaluateJacobian() == 0);

        return MapFromCsr(fixture.repca.getCsrJacobian());
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian(
          const Data& data,
          TestStatus& success,
          RealT       alpha = 1.0) const
      {
        Fixture<ScalarT> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);

        for (IdxT row = 0; row < fixture.bus.size(); ++row)
        {
          fixture.bus.setVariableIndex(row, fixture.repca.size() + row);
        }

        setAnswerKeyState(fixture.repca);
        fixture.repca.updateTime(0.0, alpha);
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= (fixture.repca.evaluateJacobian() == 0);
        success *= (fixture.repca.constructCsr() == 0);

        return MapFromCsr(fixture.repca.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
