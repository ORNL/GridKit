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
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/RepcaData.hpp>
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
    class ConverterRepcaTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterRepcaTests()  = default;
      ~ConverterRepcaTests() = default;

      // A command initialized exactly at a limit leaves the widest residual
      // tail: QPI at 708 eps.
      static constexpr RealT kBehaviorTol =
          static_cast<RealT>(1000.0) * std::numeric_limits<RealT>::epsilon();

      // The widest Enzyme-versus-dependency-tracking gap is QMEAS at the
      // qbranch column, 0.7 eps.
      static constexpr RealT kJacobianTol =
          static_cast<RealT>(10.0) * std::numeric_limits<RealT>::epsilon();

      /// Validate construction, defaults, parameters, signals, and time floors.
      TestOutcome validation()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing REPCA defaults, parameter floors, and invalid "
                         "configurations. Logged errors and warnings are expected.");

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(index(Vars::MAXIMUM)));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() > 0);

        Fixture<ScalarT> configured(makeData());
        configured.attachAllInputs();
        success *= (configured.repca.size() == static_cast<IdxT>(index(Vars::MAXIMUM)));
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

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> missing_signals(&bus, makeData());
        success *= (missing_signals.verify() > 0);

        success *= invalidParameterCase(Params::mva, 0.0);
        success *= invalidParameterCase(Params::Tfv, 0.0);
        success *= invalidParameterCase(Params::Tfv, -0.1);
        success *= invalidParameterCase(Params::dbdlow, 0.1);
        success *= invalidParameterCase(Params::dbdupper, -0.1);
        success *= invalidParameterCase(Params::emin, 0.1);
        success *= invalidParameterCase(Params::emax, -0.1);
        success *= invalidParameterCase(Params::Qmin, 1.1);
        success *= invalidParameterCase(Params::fdbd1, 0.1);
        success *= invalidParameterCase(Params::fdbd2, -0.1);
        success *= invalidParameterCase(Params::femin, 0.1);
        success *= invalidParameterCase(Params::femax, -0.1);
        success *= invalidParameterCase(Params::Pmin, 2.1);
        success *= invalidParameterCase(Params::mva, true);

        success *= invalidParameterCase(Params::Tfltr, -0.2);
        success *= invalidParameterCase(Params::Tft, -0.1);
        success *= invalidParameterCase(Params::Tp, -0.3);
        success *= invalidParameterCase(Params::Tlag, -0.4);

        const RealT nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT infinity = std::numeric_limits<RealT>::infinity();
        for (const Params parameter : {Params::mva,
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
                                       Params::Tlag})
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

        for (const Params flag : {Params::VcompFlag, Params::RefFlag, Params::Freqflag})
        {
          for (const bool value : {false, true})
          {
            auto data             = makeData();
            data.parameters[flag] = value;
            Fixture<ScalarT> model(data);
            model.attachAllInputs();
            success *= (model.repca.verify() == 0);
          }

          for (const IdxT value : {static_cast<IdxT>(0),
                                   static_cast<IdxT>(1),
                                   static_cast<IdxT>(2)})
          {
            success *= invalidParameterCase(flag, value);
          }

          for (const RealT value : {static_cast<RealT>(0.0),
                                    static_cast<RealT>(0.5),
                                    static_cast<RealT>(1.0),
                                    nan,
                                    infinity})
          {
            success *= invalidParameterCase(flag, value);
          }
        }

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> busless(nullptr, makeData());
        success *= (busless.verify() > 0);

        success *= unlinkedSignalRejected<Ext::IBRANCHR>();
        success *= unlinkedSignalRejected<Ext::IBRANCHI>();
        success *= unlinkedSignalRejected<Ext::PBRANCH>();
        success *= unlinkedSignalRejected<Ext::QBRANCH>();
        success *= unlinkedSignalRejected<Ext::FREQ>();
        success *= unlinkedSignalRejected<Ext::FREQREF>();
        success *= unlinkedSignalRejected<Ext::VREF>();
        success *= unlinkedSignalRejected<Ext::QREF>();
        success *= unlinkedSignalRejected<Ext::PPLANTREF>();

        auto floor_data                      = makeInitializationData();
        floor_data.parameters[Params::Tfltr] = 0.0;
        floor_data.parameters[Params::Tp]    = 0.0;
        floor_data.parameters[Params::Tlag]  = 0.0;

        Fixture<ScalarT> floored(floor_data);
        floored.attachAllInputs();
        setInitializationInputs(floored);
        success *= floored.initialize(0.25, 0.45);
        success *= (floored.repca.evaluateResidual() == 0);
        success *= allResidualsZero(floored.repca);

        auto* y                = floored.repca.y().getData();
        y[index(Vars::VMEAS)] -= 0.001;
        y[index(Vars::QMEAS)] -= 0.002;
        y[index(Vars::PMEAS)] -= 0.003;
        y[index(Vars::PREF)]  -= 0.004;
        floored.repca.y().setDataUpdated();
        success *= (floored.repca.evaluateResidual() == 0);
        success *= residualsMatch(floored.repca,
                                  {{Vars::VMEAS, 1.0},
                                   {Vars::QMEAS, 2.0},
                                   {Vars::PMEAS, 3.0},
                                   {Vars::PREF, 4.0}},
                                  "floored time constants");

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

        success *= stateMatches(fixture.repca,
                                {{Vars::VMEAS, 0.99200050403213},
                                 {Vars::QMEAS, 0.2},
                                 {Vars::XQPI, 0.5},
                                 {Vars::XQLAG, 0.5},
                                 {Vars::PMEAS, 0.8},
                                 {Vars::XPPI, 0.9},
                                 {Vars::PREF, 0.9},
                                 {Vars::V, 1.0},
                                 {Vars::VLDC, 0.99200050403213},
                                 {Vars::VDROOP, 1.08},
                                 {Vars::VCTRL, 0.99200050403213},
                                 {Vars::SFRZ, 1.0},
                                 {Vars::ERQ, 0.0},
                                 {Vars::ERQDB, 0.0},
                                 {Vars::ERQLIM, 0.0},
                                 {Vars::QPI, 0.5},
                                 {Vars::QEXT, 0.25},
                                 {Vars::EF, -0.00024949607977392432},
                                 {Vars::EP, 0.0},
                                 {Vars::EPLIM, 0.0},
                                 {Vars::PPI, 0.9},
                                 {Vars::PEXT, 0.45}},
                                "initialization");

        success *= scalarPreserved(fixture.input(Ext::IBRANCHR), 0.2, "preserved ibranchr");
        success *= scalarPreserved(fixture.input(Ext::IBRANCHI), -0.1, "preserved ibranchi");
        success *= scalarPreserved(fixture.input(Ext::PBRANCH), 0.4, "preserved pbranch");
        success *= scalarPreserved(fixture.input(Ext::QBRANCH), 0.1, "preserved qbranch");
        success *= scalarPreserved(fixture.input(Ext::FREQ), 0.99, "preserved freq");
        success *= scalarPreserved(fixture.qext(), 0.25, "preserved qext");
        success *= scalarPreserved(fixture.pext(), 0.45, "preserved pext");
        success *= scalarMatches(fixture.input(Ext::FREQREF), 0.99, "published freqref");
        success *= scalarMatches(fixture.input(Ext::VREF), 0.99200050403213, "published vref");
        success *= scalarMatches(fixture.input(Ext::QREF), 0.1, "published qref");
        success *= scalarMatches(fixture.input(Ext::PPLANTREF),
                                 0.39874213184871782,
                                 "published pplantref");

        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
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
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          success *= valueUnchanged(tolerances[row],
                                    absolute_tolerance,
                                    "absolute tolerance",
                                    row);
        }

        success *= monitorMatches(fixture.repca,
                                  {{0.25, 0.45, 0.99200050403213, 0.2, 0.8}},
                                  "initialization");

        success *= allResidualsZero(fixture.repca);

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
          success *= scalarPreserved(exact_commands.input(Ext::QREF), 0.1, "qref signal");
          success *= (exact_commands.repca.evaluateResidual() == 0);
          success *= allResidualsZero(exact_commands.repca);
        }

        Fixture<ScalarT> fallback(makeInitializationData(), 0.8, 0.6);
        fallback.attachRequiredInputs();
        setInitializationInputs(fallback);
        success *= fallback.initialize(0.25, 0.45);
        success *= (fallback.repca.evaluateResidual() == 0);
        success *= allResidualsZero(fallback.repca);

        Fixture<ScalarT> outputless(makeInitializationData(),
                                    0.8,
                                    0.6,
                                    100.0e6,
                                    false);
        outputless.attachAllInputs();
        setInitializationInputs(outputless);
        success *= outputless.initialize(0.25, 0.45);
        success *= (outputless.repca.evaluateResidual() == 0);
        success *= stateMatches(outputless.repca,
                                {{Vars::QEXT, 0.25},
                                 {Vars::PEXT, 0.45}},
                                "unassigned command outputs");
        success *= monitorMatches(outputless.repca,
                                  {{0.25, 0.45, 0.99200050403213, 0.2, 0.8}},
                                  "unassigned command outputs");
        success *= allResidualsZero(outputless.repca);

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
            {"line-drop/reactive-reference/disabled-frequency", true, false, false, 0.99200050403213, 0.8, 0.0},
            {"line-drop/reactive-reference/enabled-frequency", true, false, true, 0.99200050403213, 0.9, 0.45},
            {"line-drop/voltage-reference/disabled-frequency", true, true, false, 0.99200050403213, 0.8, 0.0},
            {"line-drop/voltage-reference/enabled-frequency", true, true, true, 0.99200050403213, 0.9, 0.45},
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

          success *= stateMatches(scenario.repca,
                                  {{Vars::VMEAS, test_case.voltage},
                                   {Vars::VCTRL, test_case.voltage},
                                   {Vars::PREF, test_case.pref},
                                   {Vars::PPI, test_case.pref},
                                   {Vars::PEXT, test_case.pext}},
                                  test_case.label);
          success *= scalarMatches(scenario.qext(), 0.25, test_case.label);
          success *= scalarMatches(scenario.pext(), test_case.pext, test_case.label);
          success *= allResidualsZero(scenario.repca);
        }

        return success.report(__func__);
      }

      /// Check initialization rejection, atomicity, and exact command boundaries.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible REPCA command and controller "
                         "initialization points. Logged errors are expected.");

        const auto data = makeInitializationData();

        struct RejectionCase
        {
          const char* label;
          RealT       qext;
          RealT       pext;
        };

        const RealT                        nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT                        infinity = std::numeric_limits<RealT>::infinity();
        const std::array<RejectionCase, 8> rejection_cases{{
            {"qext below Qmin", -0.4000000001, 0.45},
            {"qext above Qmax", 0.4500000001, 0.45},
            {"pext below Pmin", 0.25, -1.0e-10},
            {"pext above Pmax", 0.25, 0.6000000001},
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
        for (const Ext port : {Ext::IBRANCHR,
                               Ext::IBRANCHI,
                               Ext::PBRANCH,
                               Ext::QBRANCH,
                               Ext::FREQ})
        {
          for (const RealT value : {infinity, -infinity, nan})
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

        success *= initializationRejectedAtomically(collapsed_data,
                                                    0.2500000001,
                                                    0.45,
                                                    "qext above collapsed Q limit");
        success *= initializationRejectedAtomically(collapsed_data,
                                                    0.25,
                                                    0.4500000001,
                                                    "pext above collapsed P limit");

        auto reactive_aw_data                         = makeInitializationData();
        reactive_aw_data.parameters[Params::dbdupper] = 0.03;

        success *= initializationRejectedAtomically(reactive_aw_data,
                                                    0.25,
                                                    0.45,
                                                    "nonzero reactive-power PI antiwindup rate");
        success *= initializationRejectedAtomically(
            reactive_aw_data,
            0.25,
            0.45,
            "nonzero gated reactive-power PI state rate at the freeze transition",
            NonfiniteTarget::NONE,
            0.7,
            0.0);

        {
          Fixture<ScalarT> frozen_reactive_rate(reactive_aw_data, 0.5, 0.0);
          frozen_reactive_rate.attachAllInputs();
          setInitializationInputs(frozen_reactive_rate);
          success *= frozen_reactive_rate.initialize(0.25, 0.45);
          success *= (frozen_reactive_rate.repca.evaluateResidual() == 0);
          success *= allResidualsZero(frozen_reactive_rate.repca);
        }

        auto active_aw_data                      = makeInitializationData();
        active_aw_data.parameters[Params::femin] = 0.0;

        success *= initializationRejectedAtomically(active_aw_data,
                                                    0.25,
                                                    0.45,
                                                    "nonzero active-power PI antiwindup rate");

        auto overflow_data                    = makeInitializationData();
        overflow_data.parameters[Params::Rc]  = std::numeric_limits<RealT>::max();
        success                              *= initializationRejectedAtomically(overflow_data,
                                                    0.25,
                                                    0.45,
                                                    "nonfinite derived initialization candidate");

        // An invalid configuration is rejected before any state is written.
        {
          auto invalid_data                    = data;
          invalid_data.parameters[Params::Tfv] = 0.0;
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

        {
          Fixture<ScalarT> qmax_pmin_boundary(data, 0.8, 0.6);
          qmax_pmin_boundary.attachAllInputs();
          setInitializationInputs(qmax_pmin_boundary);
          success *= qmax_pmin_boundary.initialize(0.45, 0.0);
          success *= (qmax_pmin_boundary.repca.evaluateResidual() == 0);
          success *= stateMatches(qmax_pmin_boundary.repca,
                                  {{Vars::QPI, 0.9},
                                   {Vars::XQLAG, 0.9},
                                   {Vars::XQPI, 1.0},
                                   {Vars::QEXT, 0.45},
                                   {Vars::PREF, 0.0},
                                   {Vars::PPI, 0.0},
                                   {Vars::XPPI, -0.1},
                                   {Vars::PEXT, 0.0}},
                                  "Qmax/Pmin command boundary");
          success *= allResidualsZero(qmax_pmin_boundary.repca);
        }

        {
          Fixture<ScalarT> qmin_pmax_boundary(data, 0.8, 0.6);
          qmin_pmax_boundary.attachAllInputs();
          setInitializationInputs(qmin_pmax_boundary);
          success *= qmin_pmax_boundary.initialize(-0.4, 0.6);
          success *= (qmin_pmax_boundary.repca.evaluateResidual() == 0);
          success *= stateMatches(qmin_pmax_boundary.repca,
                                  {{Vars::QPI, -0.8},
                                   {Vars::XQLAG, -0.8},
                                   {Vars::XQPI, -0.9},
                                   {Vars::QEXT, -0.4},
                                   {Vars::PREF, 1.2},
                                   {Vars::PPI, 1.2},
                                   {Vars::XPPI, 1.3},
                                   {Vars::PEXT, 0.6}},
                                  "Qmin/Pmax command boundary");
          success *= allResidualsZero(qmin_pmax_boundary.repca);
        }

        {
          Fixture<ScalarT> collapsed_limits(collapsed_data, 0.8, 0.6);
          collapsed_limits.attachAllInputs();
          setInitializationInputs(collapsed_limits);
          success *= collapsed_limits.initialize(0.25, 0.45);
          success *= (collapsed_limits.repca.evaluateResidual() == 0);
          success *= stateMatches(collapsed_limits.repca,
                                  {{Vars::XQPI, 0.5},
                                   {Vars::XQLAG, 0.5},
                                   {Vars::QPI, 0.5},
                                   {Vars::QEXT, 0.25},
                                   {Vars::XPPI, 0.9},
                                   {Vars::PREF, 0.9},
                                   {Vars::PPI, 0.9},
                                   {Vars::PEXT, 0.45}},
                                  "collapsed Q/P limits");
          success *= allResidualsZero(collapsed_limits.repca);
        }

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

        const std::array<ExpectedResidual, index(Vars::MAXIMUM)> expected{{
            {Vars::VMEAS, "VMEAS", 0.19000000000000017},
            {Vars::QMEAS, "QMEAS", 0.77000000000000013},
            {Vars::XQPI, "XQPI", 0.018000000000000002},
            {Vars::XQLAG, "XQLAG", 0.12666666666666668},
            {Vars::PMEAS, "PMEAS", 0.92999999999999983},
            {Vars::XPPI, "XPPI", 0.082000000000000003},
            {Vars::PREF, "PREF", 0.070000000000000104},
            {Vars::V, "V", -0.02999999999999993},
            {Vars::VLDC, "VLDC", -0.015651160000000081},
            {Vars::VDROOP, "VDROOP", 0.053999999999999965},
            {Vars::VCTRL, "VCTRL", -0.030000000000000027},
            {Vars::SFRZ, "SFRZ", 0.19999999999999996},
            {Vars::ERQ, "ERQ", 2.7755575615628914e-17},
            {Vars::ERQDB, "ERQDB", -0.047111912348473055},
            {Vars::ERQLIM, "ERQLIM", 0.030000000000000044},
            {Vars::QPI, "QPI", -0.16000000000000009},
            {Vars::QEXT, "QEXT", -0.36399999999999999},
            {Vars::EF, "EF", -0.0089371400000009833},
            {Vars::EP, "EP", 0.59036181730064152},
            {Vars::EPLIM, "EPLIM", 0.049999999999999968},
            {Vars::PPI, "PPI", -0.38200000000000006},
            {Vars::PEXT, "PEXT", -0.62},
        }};

        success              *= (static_cast<size_t>(fixture.repca.getResidual().getSize()) == expected.size());
        const auto* residual  = fixture.repca.getResidual().getData();
        for (size_t row = 0; row < expected.size(); ++row)
        {
          if (index(expected[row].variable) != row)
          {
            std::cout << "REPCA residual key position " << row << " names row "
                      << expected[row].name << '\n';
            success = false;
          }
          success *= scalarMatches(residual[index(expected[row].variable)],
                                   expected[row].value,
                                   expected[row].name);
        }

        return success.report(__func__);
      }

      /// Check reactive modes, smooth limits, antiwindup, and lead-lag behavior.
      TestOutcome reactiveControl()
      {
        TestStatus success = true;

        struct FlagCase
        {
          const char* label;
          bool        voltage_compensation;
          bool        voltage_reference;
          RealT       vctrl;
          RealT       erq;
        };

        const std::array<FlagCase, 4> flag_cases{{
            {"droop/reactive-reference", false, false, 0.08, 0.30},
            {"droop/voltage-reference", false, true, 0.08, 0.10},
            {"line-drop/reactive-reference", true, false, -0.08, 0.30},
            {"line-drop/voltage-reference", true, true, -0.08, 0.10},
        }};
        for (const auto& test_case : flag_cases)
        {
          auto data                          = makeResidualData();
          data.parameters[Params::VcompFlag] = test_case.voltage_compensation;
          data.parameters[Params::RefFlag]   = test_case.voltage_reference;

          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          setAnswerKeyInputs(fixture);
          fixture.input(Ext::VREF) = 1.05;
          fixture.input(Ext::QREF) = 0.20;

          success *= fixture.prepare(0.0, 0.0);
          setState(fixture.repca,
                   {{Vars::VMEAS, 0.95},
                    {Vars::QMEAS, 0.10},
                    {Vars::VLDC, 0.92},
                    {Vars::VDROOP, 1.08},
                    {Vars::VCTRL, 1.0},
                    {Vars::ERQ, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);

          success *= residualsMatch(fixture.repca,
                                    {{Vars::VCTRL, test_case.vctrl},
                                     {Vars::ERQ, test_case.erq}},
                                    test_case.label);
        }

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);

        const std::array<DrivenCase, 3> freeze_cases{{
            {0.6, 3.7751357595539048e-11},
            {0.7, 0.5},
            {0.8, 0.99999999996224864},
        }};
        for (const auto& test_case : freeze_cases)
        {
          setState(fixture.repca, {{Vars::V, test_case.input}, {Vars::SFRZ, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca, {{Vars::SFRZ, test_case.expected}}, "freeze gate");
        }

        const std::array<DrivenCase, 5> deadband_cases{{
            {-0.05, -0.030003109594436028},
            {-0.02, -0.002888087651526948},
            {0.0, -3.104066702683552e-5},
            {0.03, 0.002888087651526948},
            {0.06, 0.030003109594436025},
        }};
        for (const auto& test_case : deadband_cases)
        {
          setState(fixture.repca, {{Vars::ERQ, test_case.input}, {Vars::ERQDB, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::ERQDB, test_case.expected}},
                                    "reactive-power deadband");
        }

        const std::array<DrivenCase, 5> error_limit_cases{{
            {-1.0, -0.7},
            {-0.7, -0.69711188674766689},
            {0.2, 0.2},
            {0.8, 0.79711188674766698},
            {1.0, 0.8},
        }};
        for (const auto& test_case : error_limit_cases)
        {
          setState(fixture.repca, {{Vars::ERQDB, test_case.input}, {Vars::ERQLIM, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::ERQLIM, test_case.expected}},
                                    "reactive-power error limit");
        }

        const std::array<DrivenCase, 5> command_limit_cases{{
            {-1.0, -0.8},
            {-0.8, -0.79711188674766698},
            {0.2, 0.2},
            {0.9, 0.89711188674766706},
            {1.1, 0.9},
        }};
        for (const auto& test_case : command_limit_cases)
        {
          setState(fixture.repca,
                   {{Vars::XQPI, test_case.input},
                    {Vars::ERQLIM, 0.0},
                    {Vars::QPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::QPI, test_case.expected}},
                                    "reactive-power command limit");
        }

        const std::array<AntiWindupCase, 8> antiwindup_cases{{
            {-0.9, -0.1, -1.1325407278661714e-11},
            {-0.9, 0.1, 0.3},
            {-0.8, -0.1, -0.15},
            {-0.8, 0.1, 0.3},
            {0.9, -0.1, -0.3},
            {0.9, 0.1, 0.15},
            {1.0, -0.1, -0.3},
            {1.0, 0.1, 1.1325407278661714e-11},
        }};
        for (const auto& test_case : antiwindup_cases)
        {
          setState(fixture.repca,
                   {{Vars::QPI, test_case.output},
                    {Vars::ERQLIM, test_case.error},
                    {Vars::SFRZ, 1.0}});
          setDerivative(fixture.repca, {{Vars::XQPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XQPI, test_case.expected}},
                                    "reactive-power antiwindup");
        }

        setState(fixture.repca,
                 {{Vars::XQLAG, 0.14},
                  {Vars::QPI, 0.27},
                  {Vars::QEXT, 0.20}});
        setDerivative(fixture.repca, {{Vars::XQLAG, -0.04}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{Vars::XQLAG, 0.12666666666666668},
                                   {Vars::QEXT, -0.36399999999999999}},
                                  "reactive-command lead-lag");

        // The command sits at Qmax with the error at emax driving further
        // out, so the antiwindup gate must block the PI state.
        {
          Fixture<DependencyTracking::Variable> blocked(makeResidualData());
          blocked.attachAllInputs();
          setAnswerKeyInputs(blocked);
          success *= blocked.prepare(0.0, 0.0);
          setState(blocked.repca,
                   {{Vars::QPI, 0.9}, {Vars::ERQLIM, 0.8}, {Vars::SFRZ, 1.0}});
          setDerivative(blocked.repca, {{Vars::XQPI, 0.0}});
          numberVariables(blocked, 1.0);
          success *= (blocked.repca.evaluateResidual() == 0);

          const JacobianRow expected{
              {index(Vars::XQPI), -1.0},
              {index(Vars::SFRZ), 1.2},
              {index(Vars::ERQLIM), 1.5},
              {index(Vars::QPI), -144.0},
          };
          success *= jacobianRowMatches(
              blocked.repca.getResidual().getData()[index(Vars::XQPI)].getDependencies(),
              expected,
              index(Vars::XQPI),
              "blocked reactive-power antiwindup",
              kBehaviorTol);
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
        success *= fixture.prepare(0.0, 0.0);

        const std::array<DrivenCase, 5> frequency_deadband_cases{{
            {-0.05, -0.040000281494001103},
            {-0.01, -0.0028777978975925616},
            {0.0, -0.00024949607977392432},
            {0.015, 0.0028777978975925616},
            {0.05, 0.035000934519396246},
        }};
        for (const auto& test_case : frequency_deadband_cases)
        {
          fixture.input(Ext::FREQ)    = 1.0;
          fixture.input(Ext::FREQREF) = 1.0 + test_case.input;
          setState(fixture.repca, {{Vars::EF, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::EF, test_case.expected}},
                                    "frequency deadband");
        }

        const std::array<DrivenCase, 3> droop_cases{{
            {-0.1, -0.099999999999842701},
            {0.0, 0.0028881132523331052},
            {0.1, 0.2000000000001573},
        }};
        fixture.input(Ext::PPLANTREF) = 0.2;
        for (const auto& test_case : droop_cases)
        {
          setState(fixture.repca,
                   {{Vars::EF, test_case.input},
                    {Vars::EP, 0.0},
                    {Vars::PMEAS, 0.4}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::EP, test_case.expected}},
                                    "frequency droop");
        }

        const std::array<DrivenCase, 5> error_limit_cases{{
            {-1.0, -0.5},
            {-0.5, -0.49711188674766688},
            {0.2, 0.2},
            {0.6, 0.59711188674766702},
            {1.0, 0.6},
        }};
        for (const auto& test_case : error_limit_cases)
        {
          setState(fixture.repca, {{Vars::EP, test_case.input}, {Vars::EPLIM, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::EPLIM, test_case.expected}},
                                    "active-power error limit");
        }

        const std::array<DrivenCase, 5> command_limit_cases{{
            {-0.2, 5.9381836780872301e-24},
            {0.0, 0.0028881132523331052},
            {0.4, 0.4},
            {1.2, 1.1971118867476669},
            {1.4, 1.2},
        }};
        for (const auto& test_case : command_limit_cases)
        {
          setState(fixture.repca,
                   {{Vars::XPPI, test_case.input},
                    {Vars::EPLIM, 0.0},
                    {Vars::PPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::PPI, test_case.expected}},
                                    "active-power command limit");
        }

        const std::array<AntiWindupCase, 8> antiwindup_cases{{
            {-0.1, -0.1, -6.7952443671970281e-12},
            {-0.1, 0.1, 0.18},
            {0.0, -0.1, -0.09},
            {0.0, 0.1, 0.18},
            {1.2, -0.1, -0.18},
            {1.2, 0.1, 0.09},
            {1.3, -0.1, -0.18},
            {1.3, 0.1, 6.7952443671970281e-12},
        }};
        for (const auto& test_case : antiwindup_cases)
        {
          setState(fixture.repca,
                   {{Vars::PPI, test_case.output},
                    {Vars::EPLIM, test_case.error}});
          setDerivative(fixture.repca, {{Vars::XPPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{Vars::XPPI, test_case.expected}},
                                    "active-power antiwindup");
        }

        setState(fixture.repca, {{Vars::PPI, 0.66}, {Vars::PREF, 0.60}});
        setDerivative(fixture.repca, {{Vars::PREF, 0.05}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{Vars::PREF, 0.070000000000000104}},
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
                      {{Vars::VMEAS, 0.11},
                       {Vars::QMEAS, -0.12},
                       {Vars::XQPI, 0.13},
                       {Vars::XQLAG, -0.14},
                       {Vars::PMEAS, 0.15},
                       {Vars::XPPI, -0.16},
                       {Vars::PREF, 0.17}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{Vars::VMEAS, 0.09000000000000018},
                                   {Vars::QMEAS, 0.8700000000000001},
                                   {Vars::XQPI, -0.082},
                                   {Vars::XQLAG, 0.22666666666666668},
                                   {Vars::PMEAS, 0.7999999999999998},
                                   {Vars::XPPI, 0.232},
                                   {Vars::PREF, -0.04999999999999989}},
                                  "explicit derivatives");

        return success.report(__func__);
      }

      /// Check every dependency-tracking Jacobian row against an independent
      /// numerical answer key, in every selector mode and at a non-unit alpha.
      TestOutcome dependencyTracking()
      {
        TestStatus success = true;

        // The fixed state places ERQ exactly at dbdupper, so the key also
        // covers the CommonMath deadband transition derivative.
        const auto data       = makeResidualData();
        const auto dependency = dependencyTrackingJacobian(data, success);

        success *= jacobianMatches(dependency,
                                   expectedJacobian(),
                                   "dependency tracking",
                                   kBehaviorTol);

        auto all_flags_off_data                          = data;
        all_flags_off_data.parameters[Params::VcompFlag] = false;
        all_flags_off_data.parameters[Params::RefFlag]   = false;
        all_flags_off_data.parameters[Params::Freqflag]  = false;
        const auto all_flags_off =
            dependencyTrackingJacobian(all_flags_off_data, success);
        success *= jacobianMatches(all_flags_off,
                                   expectedJacobianAllFlagsOff(),
                                   "all-flags-off dependency tracking",
                                   kBehaviorTol);

        const auto nonunit_alpha_dependency =
            dependencyTrackingJacobian(data, success, kNonunitAlpha);
        success *= jacobianMatches(nonunit_alpha_dependency,
                                   expectedJacobianNonunitAlpha(),
                                   "non-unit-alpha dependency tracking",
                                   kBehaviorTol);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// One rich state, every selector mode, and a non-unit alpha drive both
      /// sensitivity paths; every Enzyme CSR row must match dependency tracking.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto data = makeResidualData();

        success *= jacobianMatches(enzymeJacobian(data, success),
                                   dependencyTrackingJacobian(data, success),
                                   "Enzyme versus dependency tracking",
                                   kJacobianTol);

        auto all_flags_off_data                           = data;
        all_flags_off_data.parameters[Params::VcompFlag]  = false;
        all_flags_off_data.parameters[Params::RefFlag]    = false;
        all_flags_off_data.parameters[Params::Freqflag]   = false;
        success                                          *= jacobianMatches(
            enzymeJacobian(all_flags_off_data, success),
            dependencyTrackingJacobian(all_flags_off_data, success),
            "all-flags-off Enzyme versus dependency tracking",
            kJacobianTol);

        success *= jacobianMatches(
            enzymeJacobian(data, success, kNonunitAlpha),
            dependencyTrackingJacobian(data, success, kNonunitAlpha),
            "non-unit-alpha Enzyme versus dependency tracking",
            kJacobianTol);

        return success.report(__func__);
      }
#endif

    private:
      using Params      = PhasorDynamics::Converter::RepcaParameters;
      using Vars        = PhasorDynamics::Converter::RepcaInternalVariables;
      using Ext         = PhasorDynamics::Converter::RepcaExternalVariables;
      using Mon         = PhasorDynamics::Converter::RepcaMonitorableVariables;
      using Data        = PhasorDynamics::Converter::RepcaData<RealT, IdxT>;
      using RepcaT      = PhasorDynamics::Converter::Repca<ScalarT, IdxT>;
      using JacobianRow = DependencyTracking::Variable::DependencyMap;

      static constexpr size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(Ext variable)
      {
        return static_cast<size_t>(variable);
      }

      struct Row
      {
        constexpr Row(Vars row, RealT expected_value)
          : variable(row),
            value(expected_value)
        {
        }

        Vars  variable;
        RealT value;
      };

      struct ExpectedResidual
      {
        Vars        variable;
        const char* name;
        RealT       value;
      };

      using Rows        = std::initializer_list<Row>;
      using JacobianKey = std::array<JacobianRow, index(Vars::MAXIMUM)>;

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
        std::array<T, index(Ext::MAXIMUM)>                                   input_values_{};
        std::array<IdxT, index(Ext::MAXIMUM)>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, index(Ext::MAXIMUM)> input_nodes_{};

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
            repca.getSignals().template assignSignalNode<Vars::QEXT>(&qext_node_);
            repca.getSignals().template assignSignalNode<Vars::PEXT>(&pext_node_);
          }
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachRequiredInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = repca.size() + bus.size();
          for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = repca.getSignals();
          signals.template attachSignalNode<Ext::IBRANCHR>(&input_nodes_[index(Ext::IBRANCHR)]);
          signals.template attachSignalNode<Ext::IBRANCHI>(&input_nodes_[index(Ext::IBRANCHI)]);
          signals.template attachSignalNode<Ext::PBRANCH>(&input_nodes_[index(Ext::PBRANCH)]);
          signals.template attachSignalNode<Ext::QBRANCH>(&input_nodes_[index(Ext::QBRANCH)]);
          signals.template attachSignalNode<Ext::FREQ>(&input_nodes_[index(Ext::FREQ)]);
        }

        void attachAllInputs(RealT initial_value = 0.0)
        {
          attachRequiredInputs(initial_value);

          auto& signals = repca.getSignals();
          signals.template attachSignalNode<Ext::FREQREF>(&input_nodes_[index(Ext::FREQREF)]);
          signals.template attachSignalNode<Ext::VREF>(&input_nodes_[index(Ext::VREF)]);
          signals.template attachSignalNode<Ext::QREF>(&input_nodes_[index(Ext::QREF)]);
          signals.template attachSignalNode<Ext::PPLANTREF>(&input_nodes_[index(Ext::PPLANTREF)]);
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

        PhasorDynamics::Bus<T, IdxT>              bus;
        PhasorDynamics::Converter::Repca<T, IdxT> repca;
      };

      static constexpr RealT kStateVr      = 0.9;
      static constexpr RealT kStateVi      = 0.4;
      static constexpr RealT kNonunitAlpha = 2.5;

      static constexpr size_t kBusVrColumn        = index(Vars::MAXIMUM);
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
        data.parameters[Params::Tfv]      = 1.5;
        data.parameters[Params::Tp]       = 0.4;
        data.parameters[Params::fdbd1]    = -0.01;
        data.parameters[Params::fdbd2]    = 0.015;
        data.parameters[Params::Ddn]      = 2.0;
        data.parameters[Params::Dup]      = 1.0;
        data.parameters[Params::femax]    = 0.6;
        data.parameters[Params::femin]    = -0.5;
        data.parameters[Params::Kpg]      = 1.7;
        data.parameters[Params::Kig]      = 1.8;
        data.parameters[Params::Pmax]     = 1.2;
        data.parameters[Params::Tlag]     = 0.5;
        return data;
      }

      Data makeInitializationData() const
      {
        auto data                         = makeResidualData();
        data.parameters[Params::dbdupper] = 0.02;
        data.parameters[Params::emin]     = -0.8;
        data.parameters[Params::femin]    = -0.6;
        return data;
      }

      template <typename T>
      void setInitializationInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::IBRANCHR) = static_cast<T>(0.2);
        fixture.input(Ext::IBRANCHI) = static_cast<T>(-0.1);
        fixture.input(Ext::PBRANCH)  = static_cast<T>(0.4);
        fixture.input(Ext::QBRANCH)  = static_cast<T>(0.1);
        fixture.input(Ext::FREQ)     = static_cast<T>(0.99);
      }

      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::IBRANCHR)  = static_cast<T>(0.08);
        fixture.input(Ext::IBRANCHI)  = static_cast<T>(-0.02);
        fixture.input(Ext::PBRANCH)   = static_cast<T>(0.41);
        fixture.input(Ext::QBRANCH)   = static_cast<T>(0.13);
        fixture.input(Ext::FREQ)      = static_cast<T>(0.99);
        fixture.input(Ext::FREQREF)   = static_cast<T>(1.0);
        fixture.input(Ext::VREF)      = static_cast<T>(1.01);
        fixture.input(Ext::QREF)      = static_cast<T>(0.12);
        fixture.input(Ext::PPLANTREF) = static_cast<T>(0.55);
      }

      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Converter::Repca<T, IdxT>& repca) const
      {
        setState(repca,
                 {{Vars::VMEAS, 0.98},
                  {Vars::QMEAS, 0.11},
                  {Vars::XQPI, 0.07},
                  {Vars::XQLAG, 0.14},
                  {Vars::PMEAS, 0.44},
                  {Vars::XPPI, 0.21},
                  {Vars::PREF, 0.60},
                  {Vars::V, 1.0},
                  {Vars::VLDC, 0.99},
                  {Vars::VDROOP, 1.05},
                  {Vars::VCTRL, 1.02},
                  {Vars::SFRZ, 0.8},
                  {Vars::ERQ, 0.03},
                  {Vars::ERQDB, 0.05},
                  {Vars::ERQLIM, 0.02},
                  {Vars::QPI, 0.27},
                  {Vars::QEXT, 0.20},
                  {Vars::EF, 0.01},
                  {Vars::EP, 0.09},
                  {Vars::EPLIM, 0.04},
                  {Vars::PPI, 0.66},
                  {Vars::PEXT, 0.61}});
        setDerivative(repca,
                      {{Vars::VMEAS, 0.01},
                       {Vars::QMEAS, -0.02},
                       {Vars::XQPI, 0.03},
                       {Vars::XQLAG, -0.04},
                       {Vars::PMEAS, 0.02},
                       {Vars::XPPI, -0.01},
                       {Vars::PREF, 0.05}});
      }

      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), 0.9, 0.4);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), 0.9, 0.4);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        implicit_defaults.input(Ext::PBRANCH) = 0.2;
        implicit_defaults.input(Ext::QBRANCH) = 0.1;
        implicit_defaults.input(Ext::FREQ)    = 1.0;
        explicit_defaults.input(Ext::PBRANCH) = 0.2;
        explicit_defaults.input(Ext::QBRANCH) = 0.1;
        explicit_defaults.input(Ext::FREQ)    = 1.0;

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
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
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
        fixture.repca.getSignals().template attachSignalNode<variable>(&unlinked_node);
        return fixture.repca.verify() > 0;
      }

      /// Fill state and derivative with a recognizable ramp, restoring the
      /// aliased commands, so any write by a rejected initialization shows.
      void poisonState(Fixture<ScalarT>& fixture, RealT qext, RealT pext) const
      {
        auto* y  = fixture.repca.y().getData();
        auto* yp = fixture.repca.yp().getData();
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
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
                                            Ext             poisoned_port = Ext::FREQ,
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

        const auto                             y_before   = copyVector(fixture.repca.y());
        const auto                             yp_before  = copyVector(fixture.repca.yp());
        const auto                             bus_before = copyVector(fixture.bus.y());
        std::array<RealT, index(Ext::MAXIMUM)> inputs_before{};
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
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
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
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
      void setState(PhasorDynamics::Converter::Repca<T, IdxT>& repca, Rows rows) const
      {
        auto* y = repca.y().getData();
        for (const auto& [variable, value] : rows)
        {
          y[index(variable)] = static_cast<T>(value);
        }
        repca.y().setDataUpdated();
      }

      template <typename T>
      void setDerivative(PhasorDynamics::Converter::Repca<T, IdxT>& repca, Rows rows) const
      {
        auto* yp = repca.yp().getData();
        for (const auto& [variable, value] : rows)
        {
          yp[index(variable)] = static_cast<T>(value);
        }
        repca.yp().setDataUpdated();
      }

      static const char* variableName(Vars variable)
      {
        static constexpr std::array<const char*, index(Vars::MAXIMUM)> names{{
            "VMEAS",
            "QMEAS",
            "XQPI",
            "XQLAG",
            "PMEAS",
            "XPPI",
            "PREF",
            "V",
            "VLDC",
            "VDROOP",
            "VCTRL",
            "SFRZ",
            "ERQ",
            "ERQDB",
            "ERQLIM",
            "QPI",
            "QEXT",
            "EF",
            "EP",
            "EPLIM",
            "PPI",
            "PEXT",
        }};
        return names[index(variable)];
      }

      static bool variableMatches(RealT       actual,
                                  RealT       expected,
                                  const char* what,
                                  Vars        variable,
                                  const char* context,
                                  RealT       tolerance = kBehaviorTol)
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
                             RealT       tolerance = kBehaviorTol)
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
                         RealT       tolerance = kBehaviorTol) const
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

      template <typename VectorT>
      bool rowsMatch(const VectorT& vector,
                     const Row*     rows,
                     size_t         count,
                     const char*    what,
                     const char*    context) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t i = 0; i < count; ++i)
        {
          const auto& [variable, expected] = rows[i];
          const size_t row                 = index(variable);

          if (!variableMatches(static_cast<RealT>(values[row]),
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

      bool residualsMatch(const RepcaT& repca, Rows rows, const char* context = "") const
      {
        return rowsMatch(repca.getResidual(),
                         rows.begin(),
                         rows.size(),
                         "residual",
                         context);
      }

      bool stateMatches(const RepcaT& repca, Rows rows, const char* context = "") const
      {
        return rowsMatch(repca.y(), rows.begin(), rows.size(), "state", context);
      }

      bool allResidualsZero(const RepcaT& repca) const
      {
        bool        success = true;
        const auto* f       = repca.getResidual().getData();
        const auto* yp      = repca.yp().getData();
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          const auto variable = static_cast<Vars>(row);
          if (!variableMatches(f[row],
                               0.0,
                               "residual",
                               variable,
                               "at rest"))
          {
            success = false;
          }
          if (!variableMatches(yp[row],
                               0.0,
                               "derivative",
                               variable,
                               "at rest"))
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

      void noteExpectedLogs(const char* message) const
      {
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << '\n';
        Log::setVerbosity(previous_verbosity);
      }

      JacobianKey expectedJacobian() const
      {
        return {{
            {{index(Vars::VMEAS), -6.0}, {index(Vars::VCTRL), 5.0}},
            {{index(Vars::QMEAS), -6.0}, {externalColumn(index(Ext::QBRANCH)), 10.0}},
            {{index(Vars::XQPI), -1.0}, {index(Vars::SFRZ), 0.06}, {index(Vars::ERQLIM), 2.4}},
            {{index(Vars::XQLAG), -1.6666666666666667}, {index(Vars::QPI), 0.6666666666666666}},
            {{index(Vars::PMEAS), -3.5}, {externalColumn(index(Ext::PBRANCH)), 5.0}},
            {{index(Vars::XPPI), -1.0}, {index(Vars::EPLIM), 1.8}},
            {{index(Vars::PREF), -3.0}, {index(Vars::PPI), 2.0}},
            {{index(Vars::V), -2.0}, {kBusVrColumn, 1.8}, {kBusViColumn, 0.8}},
            {{index(Vars::VLDC), -1.98},
             {kBusVrColumn, 1.7956},
             {kBusViColumn, 0.796},
             {externalColumn(index(Ext::IBRANCHR)), -0.059792},
             {externalColumn(index(Ext::IBRANCHI)), 0.037948}},
            {{index(Vars::V), 1.0}, {index(Vars::VDROOP), -1.0}, {externalColumn(index(Ext::QBRANCH)), 0.8}},
            {{index(Vars::VLDC), 1.0}, {index(Vars::VCTRL), -1.0}},
            {{index(Vars::SFRZ), -1.0}},
            {{index(Vars::VMEAS), -1.0}, {index(Vars::ERQ), -1.0}, {externalColumn(index(Ext::VREF)), 1.0}},
            {{index(Vars::ERQ), 0.50000614417460221}, {index(Vars::ERQDB), -1.0}},
            {{index(Vars::ERQDB), 1.0}, {index(Vars::ERQLIM), -1.0}},
            {{index(Vars::XQPI), 1.0}, {index(Vars::ERQLIM), 2.0}, {index(Vars::QPI), -1.0}},
            {{index(Vars::XQLAG), 1.3}, {index(Vars::QPI), 0.2}, {index(Vars::QEXT), -3.0}},
            {{index(Vars::EF), -1.0},
             {externalColumn(index(Ext::FREQ)), -0.23963778765414229},
             {externalColumn(index(Ext::FREQREF)), 0.23963778765414229}},
            {{index(Vars::PMEAS), -1.0},
             {index(Vars::EF), 1.9168273035060777},
             {index(Vars::EP), -1.0},
             {externalColumn(index(Ext::PPLANTREF)), 2.0}},
            {{index(Vars::EP), 1.0}, {index(Vars::EPLIM), -1.0}},
            {{index(Vars::XPPI), 1.0}, {index(Vars::EPLIM), 1.7}, {index(Vars::PPI), -1.0}},
            {{index(Vars::PREF), 1.0}, {index(Vars::PEXT), -2.0}},
        }};
      }

      JacobianKey expectedJacobianAllFlagsOff() const
      {
        auto expected                = expectedJacobian();
        expected[index(Vars::VCTRL)] = {
            {index(Vars::VDROOP), 1.0},
            {index(Vars::VCTRL), -1.0},
        };
        expected[index(Vars::ERQ)] = {
            {index(Vars::QMEAS), -1.0},
            {index(Vars::ERQ), -1.0},
            {externalColumn(index(Ext::QREF)), 2.0},
        };
        expected[index(Vars::PEXT)] = {
            {index(Vars::PEXT), -2.0},
        };
        return expected;
      }

      JacobianKey expectedJacobianNonunitAlpha() const
      {
        auto expected                                    = expectedJacobian();
        expected[index(Vars::VMEAS)][index(Vars::VMEAS)] = -7.5;
        expected[index(Vars::QMEAS)][index(Vars::QMEAS)] = -7.5;
        expected[index(Vars::XQPI)][index(Vars::XQPI)]   = -2.5;
        expected[index(Vars::XQLAG)][index(Vars::XQLAG)] = -3.1666666666666665;
        expected[index(Vars::PMEAS)][index(Vars::PMEAS)] = -5.0;
        expected[index(Vars::XPPI)][index(Vars::XPPI)]   = -2.5;
        expected[index(Vars::PREF)][index(Vars::PREF)]   = -4.5;
        return expected;
      }

      static RealT mapValue(const JacobianRow& row, size_t column)
      {
        const auto entry = row.find(column);
        if (entry == row.end())
        {
          return 0.0;
        }
        return entry->second;
      }

      bool jacobianRowMatches(const JacobianRow& actual,
                              const JacobianRow& expected,
                              size_t             row,
                              const char*        source,
                              RealT              tolerance) const
      {
        constexpr size_t column_count = externalColumn(index(Ext::MAXIMUM));
        bool             success      = true;
        for (size_t column = 0; column < column_count; ++column)
        {
          const RealT actual_value   = mapValue(actual, column);
          const RealT expected_value = mapValue(expected, column);
          if (!isEqual(actual_value, expected_value, tolerance))
          {
            std::cout << "REPCA " << source << " Jacobian row " << row
                      << " column " << column << " mismatch: "
                      << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual_value << " != "
                      << expected_value << '\n';
            success = false;
          }
        }
        for (const auto& [column, value] : actual)
        {
          if (column >= column_count && !isEqual(value, 0.0, tolerance))
          {
            std::cout << "REPCA " << source << " Jacobian row " << row
                      << " has unexpected column " << column << '\n';
            success = false;
          }
        }
        return success;
      }

      template <typename ExpectedT>
      bool jacobianMatches(const std::vector<JacobianRow>& actual,
                           const ExpectedT&                expected,
                           const char*                     source,
                           RealT                           tolerance) const
      {
        if (actual.size() != expected.size())
        {
          std::cout << "REPCA " << source << " Jacobian row-count mismatch\n";
          return false;
        }

        bool success = true;
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          if (!jacobianRowMatches(actual[row], expected[row], row, source, tolerance))
          {
            success = false;
          }
        }
        return success;
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture,
                           RealT                                  alpha) const
      {
        auto* y     = fixture.repca.y().getData();
        auto* yp    = fixture.repca.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          y[row].setVariableNumber(row);
          yp[row].setVariableNumber(row);
          yp[row].scaleDependencies(alpha);
        }
        for (size_t row = 0; row < static_cast<size_t>(fixture.bus.size()); ++row)
        {
          bus_y[row].setVariableNumber(kBusVrColumn + row);
        }
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
          const auto variable = static_cast<Ext>(port);
          fixture.input(variable).setVariableNumber(fixture.inputIndex(variable));
        }

        fixture.repca.y().setDataUpdated();
        fixture.repca.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<JacobianRow> dependencyTrackingJacobian(const Data& data,
                                                          TestStatus& success,
                                                          RealT       alpha = 1.0) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        numberVariables(fixture, alpha);
        success *= (fixture.repca.evaluateResidual() == 0);

        std::vector<JacobianRow> rows(index(Vars::MAXIMUM));
        const auto*              f = fixture.repca.getResidual().getData();
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          rows[row] = f[row].getDependencies();
        }
        return rows;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      std::vector<JacobianRow> enzymeJacobian(const Data& data,
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
