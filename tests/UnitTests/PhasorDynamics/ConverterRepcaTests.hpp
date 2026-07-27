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

      // Exact-limit limiter seeds leave residual tails below 2e-13. This margin
      // absorbs them while resolving the approximately 6.8e-12 boundary response.
      static constexpr RealT kBehaviorTol = 1.0e-12;

      // Enzyme traverses the smooth expressions through a separate AD path;
      // its double-precision derivatives agree with the fixed key to O(1e-10).
      static constexpr RealT kJacobianTol = 1.0e-9;

      /// Validate construction, defaults, parameters, signals, and time floors.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (empty.getMonitor() == nullptr);

        Fixture<ScalarT> configured(makeData());
        configured.attachAllInputs();
        success *= (configured.repca.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (configured.repca.getMonitor() != nullptr);
        success *= (configured.repca.verify() == 0);

        noteExpectedLogs("Testing REPCA defaults, parameter floors, and invalid "
                         "configurations. Logged errors and warnings are expected.");

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

        for (const Params flag : {Params::VcompFlag, Params::RefFlag, Params::Freqflag})
        {
          for (const auto value : {static_cast<IdxT>(0), static_cast<IdxT>(1)})
          {
            auto data             = makeData();
            data.parameters[flag] = value;
            Fixture<ScalarT> model(data);
            model.attachAllInputs();
            success *= (model.repca.verify() == 0);
          }

          for (const auto value : {static_cast<RealT>(0.0), static_cast<RealT>(1.0)})
          {
            auto data             = makeData();
            data.parameters[flag] = value;
            Fixture<ScalarT> model(data);
            model.attachAllInputs();
            success *= (model.repca.verify() == 0);
          }

          success *= invalidParameterCase(flag, static_cast<IdxT>(2));
          success *= invalidParameterCase(flag, static_cast<RealT>(0.5));
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
        floor_data.parameters[Params::Tfltr] = -0.2;
        floor_data.parameters[Params::Tp]    = 0.0;
        floor_data.parameters[Params::Tlag]  = -0.4;

        Fixture<ScalarT> floored(floor_data);
        floored.attachAllInputs();
        setInitializationInputs(floored);
        success *= floored.initialize(0.25, 0.45);
        success *= (floored.repca.evaluateResidual() == 0);
        success *= allResidualsZero(floored.repca);

        auto* y      = floored.repca.y().getData();
        y[I::VMEAS] -= 0.001;
        y[I::QMEAS] -= 0.002;
        y[I::PMEAS] -= 0.003;
        y[I::PREF]  -= 0.004;
        floored.repca.y().setDataUpdated();
        success *= (floored.repca.evaluateResidual() == 0);
        success *= residualsMatch(floored.repca,
                                  {{I::VMEAS, 1.0},
                                   {I::QMEAS, 2.0},
                                   {I::PMEAS, 3.0},
                                   {I::PREF, 4.0}},
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
                                {{I::VMEAS, 0.99200050403213},
                                 {I::QMEAS, 0.2},
                                 {I::XQPI, 0.5},
                                 {I::XQLAG, 0.5},
                                 {I::PMEAS, 0.8},
                                 {I::XPPI, 0.9},
                                 {I::PREF, 0.9},
                                 {I::V, 1.0},
                                 {I::VLDC, 0.99200050403213},
                                 {I::VDROOP, 1.08},
                                 {I::VCTRL, 0.99200050403213},
                                 {I::SFRZ, 1.0},
                                 {I::ERQ, 0.0},
                                 {I::ERQDB, 0.0},
                                 {I::ERQLIM, 0.0},
                                 {I::QPI, 0.5},
                                 {I::QEXT, 0.25},
                                 {I::EF, -0.00024949607977392432},
                                 {I::EP, 0.0},
                                 {I::EPLIM, 0.0},
                                 {I::PPI, 0.9},
                                 {I::PEXT, 0.45}},
                                "initialization");

        success *= scalarMatches(fixture.input(E::IBRANCHR), 0.2, "preserved ibranchr");
        success *= scalarMatches(fixture.input(E::IBRANCHI), -0.1, "preserved ibranchi");
        success *= scalarMatches(fixture.input(E::PBRANCH), 0.4, "preserved pbranch");
        success *= scalarMatches(fixture.input(E::QBRANCH), 0.1, "preserved qbranch");
        success *= scalarMatches(fixture.input(E::FREQ), 0.99, "preserved freq");
        success *= scalarMatches(fixture.qext(), 0.25, "preserved qext");
        success *= scalarMatches(fixture.pext(), 0.45, "preserved pext");
        success *= scalarMatches(fixture.input(E::FREQREF), 0.99, "published freqref");
        success *= scalarMatches(fixture.input(E::VREF), 0.99200050403213, "published vref");
        success *= scalarMatches(fixture.input(E::QREF), 0.1, "published qref");
        success *= scalarMatches(fixture.input(E::PPLANTREF),
                                 0.39874213184871782,
                                 "published pplantref");

        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          const bool expected = row <= I::PREF;
          if (fixture.repca.tag()[row] != expected)
          {
            std::cout << "REPCA differentiability tag " << row << " mismatch\n";
            success = false;
          }
        }

        constexpr RealT absolute_tolerance = 2.5e-7;

        success *= (fixture.repca.setAbsoluteTolerance(absolute_tolerance) == 0);

        const auto* tolerances = fixture.repca.absoluteTolerance().getData();
        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          success *= valueUnchanged(tolerances[row],
                                    absolute_tolerance,
                                    "absolute tolerance",
                                    row);
        }

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.repca.getMonitor());
        std::stringstream monitor_output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, monitor_output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string monitor_header;
        std::string monitor_values;
        std::getline(monitor_output, monitor_header);
        std::getline(monitor_output, monitor_values);
        success *= (monitor_header == "t,Repca_repca_test_qext,Repca_repca_test_pext,"
                                      "Repca_repca_test_vmeas,Repca_repca_test_qmeas,"
                                      "Repca_repca_test_pmeas");

        const auto monitored = Tokenizer<RealT>(monitor_values, ',')();
        if (monitored.size() == 6)
        {
          success *= scalarMatches(monitored[1], 0.25, "monitored qext");
          success *= scalarMatches(monitored[2], 0.45, "monitored pext");
          success *= scalarMatches(monitored[3], 0.99200050403213, "monitored vmeas");
          success *= scalarMatches(monitored[4], 0.2, "monitored qmeas");
          success *= scalarMatches(monitored[5], 0.8, "monitored pmeas");
        }
        else
        {
          std::cout << "REPCA monitor emitted " << monitored.size()
                    << " values instead of 6\n";
          success = false;
        }

        success *= allResidualsZero(fixture.repca);

        {
          auto exact_data                    = makeInitializationData();
          exact_data.parameters[Params::mva] = 73.0;

          // This non-binary base ratio exposes any component-base round trip.
          Fixture<ScalarT> exact_commands(exact_data, 0.8, 0.6);
          exact_commands.attachAllInputs();
          setInitializationInputs(exact_commands);
          success *= exact_commands.initialize(0.25, 0.45);
          success *= valueUnchanged(exact_commands.qext(), 0.25, "qext signal", I::QEXT);
          success *= valueUnchanged(exact_commands.pext(), 0.45, "pext signal", I::PEXT);
          success *= valueUnchanged(exact_commands.input(E::QREF),
                                    0.1,
                                    "qref signal",
                                    E::QREF);
          success *= (exact_commands.repca.evaluateResidual() == 0);
          success *= allResidualsZero(exact_commands.repca);
        }

        Fixture<ScalarT> fallback(makeInitializationData(), 0.8, 0.6);
        fallback.attachRequiredInputs();
        setInitializationInputs(fallback);
        success *= fallback.initialize(0.25, 0.45);
        success *= (fallback.repca.evaluateResidual() == 0);
        success *= allResidualsZero(fallback.repca);

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
                                  {{I::VMEAS, test_case.voltage},
                                   {I::VCTRL, test_case.voltage},
                                   {I::PREF, test_case.pref},
                                   {I::PPI, test_case.pref},
                                   {I::PEXT, test_case.pext}},
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

        const RealT                        infinity = std::numeric_limits<RealT>::infinity();
        const std::array<RejectionCase, 6> rejection_cases{{
            {"qext below Qmin", -0.4000000001, 0.45},
            {"qext above Qmax", 0.4500000001, 0.45},
            {"pext below Pmin", 0.25, -1.0e-10},
            {"pext above Pmax", 0.25, 0.6000000001},
            {"nonfinite qext", infinity, 0.45},
            {"nonfinite pext", 0.25, infinity},
        }};
        for (const auto& test_case : rejection_cases)
        {
          success *= initializationRejectedAtomically(data,
                                                      test_case.qext,
                                                      test_case.pext,
                                                      test_case.label);
        }
        success *= initializationRejectedAtomically(data,
                                                    0.25,
                                                    0.45,
                                                    "nonfinite required signal",
                                                    NonfiniteTarget::FREQUENCY);
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

        auto active_aw_data                      = makeInitializationData();
        active_aw_data.parameters[Params::femin] = 0.0;

        success *= initializationRejectedAtomically(active_aw_data,
                                                    0.25,
                                                    0.45,
                                                    "nonzero active-power PI antiwindup rate");

        {
          Fixture<ScalarT> qmax_pmin_boundary(data, 0.8, 0.6);
          qmax_pmin_boundary.attachAllInputs();
          setInitializationInputs(qmax_pmin_boundary);
          success *= qmax_pmin_boundary.initialize(0.45, 0.0);
          success *= (qmax_pmin_boundary.repca.evaluateResidual() == 0);
          success *= stateMatches(qmax_pmin_boundary.repca,
                                  {{I::QPI, 0.9},
                                   {I::XQLAG, 0.9},
                                   {I::XQPI, 1.0},
                                   {I::QEXT, 0.45},
                                   {I::PREF, 0.0},
                                   {I::PPI, 0.0},
                                   {I::XPPI, -0.1},
                                   {I::PEXT, 0.0}},
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
                                  {{I::QPI, -0.8},
                                   {I::XQLAG, -0.8},
                                   {I::XQPI, -0.9},
                                   {I::QEXT, -0.4},
                                   {I::PREF, 1.2},
                                   {I::PPI, 1.2},
                                   {I::XPPI, 1.3},
                                   {I::PEXT, 0.6}},
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
                                  {{I::XQPI, 0.5},
                                   {I::XQLAG, 0.5},
                                   {I::QPI, 0.5},
                                   {I::QEXT, 0.25},
                                   {I::XPPI, 0.9},
                                   {I::PREF, 0.9},
                                   {I::PPI, 0.9},
                                   {I::PEXT, 0.45}},
                                  "collapsed Q/P limits");
          success *= allResidualsZero(collapsed_limits.repca);
        }

        return success.report(__func__);
      }

      /// Check every residual row against an independent numerical answer key.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeDynamicData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        success *= (fixture.repca.evaluateResidual() == 0);

        const std::array<Row, I::MAXIMUM> expected{{
            {I::VMEAS, 0.19000000000000017},
            {I::QMEAS, 0.77000000000000013},
            {I::XQPI, 0.018000000000000002},
            {I::XQLAG, 0.12666666666666668},
            {I::PMEAS, 0.92999999999999983},
            {I::XPPI, 0.082000000000000003},
            {I::PREF, 0.070000000000000104},
            {I::V, -0.029999999999999916},
            {I::VLDC, -0.015651160000000081},
            {I::VDROOP, 0.053999999999999965},
            {I::VCTRL, -0.030000000000000027},
            {I::SFRZ, 0.19999999999999996},
            {I::ERQ, 2.7755575615628914e-17},
            {I::ERQDB, -0.017111912348473052},
            {I::ERQLIM, 1.7347234759768071e-17},
            {I::QPI, -0.16000000000000009},
            {I::QEXT, -0.36399999999999999},
            {I::EF, -0.0089371400000009833},
            {I::EP, 0.64036181730064157},
            {I::EPLIM, 3.4694469519536142e-17},
            {I::PPI, -0.38200000000000001},
            {I::PEXT, -0.62},
        }};

        success *= (static_cast<size_t>(fixture.repca.getResidual().getSize()) == expected.size());
        for (size_t row = 0; row < expected.size(); ++row)
        {
          if (expected[row].first != row)
          {
            std::cout << "REPCA residual key position " << row << " names row "
                      << expected[row].first << '\n';
            success = false;
          }
        }
        success *= residualsMatch(fixture.repca, expected);

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
          auto data                          = makeDynamicData();
          data.parameters[Params::VcompFlag] = test_case.voltage_compensation;
          data.parameters[Params::RefFlag]   = test_case.voltage_reference;

          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          setAnswerKeyInputs(fixture);
          fixture.input(E::VREF) = 1.05;
          fixture.input(E::QREF) = 0.20;

          success *= fixture.prepare(0.0, 0.0);
          setState(fixture.repca,
                   {{I::VMEAS, 0.95},
                    {I::QMEAS, 0.10},
                    {I::VLDC, 0.92},
                    {I::VDROOP, 1.08},
                    {I::VCTRL, 1.0},
                    {I::ERQ, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);

          success *= residualsMatch(fixture.repca,
                                    {{I::VCTRL, test_case.vctrl},
                                     {I::ERQ, test_case.erq}},
                                    test_case.label);
        }

        Fixture<ScalarT> fixture(makeDynamicData());
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
          setState(fixture.repca, {{I::V, test_case.input}, {I::SFRZ, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca, {{I::SFRZ, test_case.expected}}, "freeze gate");
        }

        const std::array<DrivenCase, 4> deadband_cases{{
            {-0.05, -0.030003109594436028},
            {0.0, -3.104066702683552e-5},
            {0.03, 0.002888087651526948},
            {0.06, 0.030003109594436025},
        }};
        for (const auto& test_case : deadband_cases)
        {
          setState(fixture.repca, {{I::ERQ, test_case.input}, {I::ERQDB, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::ERQDB, test_case.expected}},
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
          setState(fixture.repca, {{I::ERQDB, test_case.input}, {I::ERQLIM, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::ERQLIM, test_case.expected}},
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
                   {{I::XQPI, test_case.input},
                    {I::ERQLIM, 0.0},
                    {I::QPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::QPI, test_case.expected}},
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
                   {{I::QPI, test_case.output},
                    {I::ERQLIM, test_case.error},
                    {I::SFRZ, 1.0}});
          setDerivative(fixture.repca, {{I::XQPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::XQPI, test_case.expected}},
                                    "reactive-power antiwindup");
        }

        setState(fixture.repca,
                 {{I::XQLAG, 0.14},
                  {I::QPI, 0.27},
                  {I::QEXT, 0.20}});
        setDerivative(fixture.repca, {{I::XQLAG, -0.04}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{I::XQLAG, 0.12666666666666668},
                                   {I::QEXT, -0.36399999999999999}},
                                  "reactive-command lead-lag");

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
          auto data                         = makeDynamicData();
          data.parameters[Params::Freqflag] = test_case.frequency_control;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          setAnswerKeyInputs(fixture);
          success *= fixture.prepare(0.0, 0.0);
          setState(fixture.repca, {{I::PREF, 0.8}, {I::PEXT, 0.3}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::PEXT, test_case.pext}},
                                    test_case.label);
        }

        Fixture<ScalarT> fixture(makeDynamicData());
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);

        const std::array<DrivenCase, 4> frequency_deadband_cases{{
            {-0.05, -0.040000281494001103},
            {0.0, -0.00024949607977392432},
            {0.015, 0.0028777978975925616},
            {0.05, 0.035000934519396246},
        }};
        for (const auto& test_case : frequency_deadband_cases)
        {
          fixture.input(E::FREQ)    = 1.0;
          fixture.input(E::FREQREF) = 1.0 + test_case.input;
          setState(fixture.repca, {{I::EF, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::EF, test_case.expected}},
                                    "frequency deadband");
        }

        const std::array<DrivenCase, 3> droop_cases{{
            {-0.1, -0.099999999999842701},
            {0.0, 0.0028881132523331052},
            {0.1, 0.2000000000001573},
        }};
        fixture.input(E::PPLANTREF) = 0.2;
        for (const auto& test_case : droop_cases)
        {
          setState(fixture.repca,
                   {{I::EF, test_case.input},
                    {I::EP, 0.0},
                    {I::PMEAS, 0.4}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::EP, test_case.expected}},
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
          setState(fixture.repca, {{I::EP, test_case.input}, {I::EPLIM, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::EPLIM, test_case.expected}},
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
                   {{I::XPPI, test_case.input},
                    {I::EPLIM, 0.0},
                    {I::PPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::PPI, test_case.expected}},
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
                   {{I::PPI, test_case.output},
                    {I::EPLIM, test_case.error}});
          setDerivative(fixture.repca, {{I::XPPI, 0.0}});
          success *= (fixture.repca.evaluateResidual() == 0);
          success *= residualsMatch(fixture.repca,
                                    {{I::XPPI, test_case.expected}},
                                    "active-power antiwindup");
        }

        setState(fixture.repca, {{I::PPI, 0.66}, {I::PREF, 0.60}});
        setDerivative(fixture.repca, {{I::PREF, 0.05}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{I::PREF, 0.070000000000000104}},
                                  "active-power command lag");

        return success.report(__func__);
      }

      /// Check every differential residual with nonzero explicit derivatives.
      TestOutcome derivatives()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeDynamicData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        setDerivative(fixture.repca,
                      {{I::VMEAS, 0.11},
                       {I::QMEAS, -0.12},
                       {I::XQPI, 0.13},
                       {I::XQLAG, -0.14},
                       {I::PMEAS, 0.15},
                       {I::XPPI, -0.16},
                       {I::PREF, 0.17}});
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= residualsMatch(fixture.repca,
                                  {{I::VMEAS, 0.09000000000000018},
                                   {I::QMEAS, 0.8700000000000001},
                                   {I::XQPI, -0.082},
                                   {I::XQLAG, 0.22666666666666668},
                                   {I::PMEAS, 0.7999999999999998},
                                   {I::XPPI, 0.232},
                                   {I::PREF, -0.04999999999999989}},
                                  "explicit derivatives");

        return success.report(__func__);
      }

      /// Check every Jacobian row against an independent numerical answer key.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto data       = makeDynamicData();
        const auto expected   = expectedJacobian();
        const auto dependency = dependencyTrackingJacobian(data, success);

        success *= jacobianMatches(dependency,
                                   expected,
                                   "dependency tracking",
                                   kBehaviorTol);

#ifdef GRIDKIT_ENABLE_ENZYME
        const auto enzyme = enzymeJacobian(data, success);

        success *= jacobianMatches(enzyme, expected, "Enzyme", kJacobianTol);
#endif

        return success.report(__func__);
      }

    private:
      using Params      = PhasorDynamics::Converter::RepcaParameters;
      using Vars        = PhasorDynamics::Converter::RepcaInternalVariables;
      using Ext         = PhasorDynamics::Converter::RepcaExternalVariables;
      using Mon         = PhasorDynamics::Converter::RepcaMonitorableVariables;
      using Data        = PhasorDynamics::Converter::RepcaData<RealT, IdxT>;
      using I           = PhasorDynamics::Converter::RepcaIdx;
      using E           = PhasorDynamics::Converter::RepcaExt;
      using RepcaT      = PhasorDynamics::Converter::Repca<ScalarT, IdxT>;
      using Row         = std::pair<size_t, RealT>;
      using Rows        = std::initializer_list<Row>;
      using JacobianRow = DependencyTracking::Variable::DependencyMap;

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
        FREQUENCY,
        BUS_VOLTAGE
      };

      /// Owns the regulated bus, REPCA, assigned command nodes, and attached
      /// input nodes. Signal storage precedes the model so every referenced node
      /// outlives REPCA; copying would invalidate the model and node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, E::MAXIMUM>                                   input_values_{};
        std::array<IdxT, E::MAXIMUM>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, E::MAXIMUM> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> qext_node_;
        PhasorDynamics::SignalNode<T, IdxT> pext_node_;

      public:
        explicit Fixture(const Data& data,
                         RealT       vr             = 1.0,
                         RealT       vi             = 0.0,
                         RealT       system_va_base = 100.0e6)
          : bus(static_cast<T>(vr), static_cast<T>(vi)),
            repca(&bus, data)
        {
          repca.setSystemBase(60.0, system_va_base);
          repca.getSignals().template assignSignalNode<Vars::QEXT>(&qext_node_);
          repca.getSignals().template assignSignalNode<Vars::PEXT>(&pext_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachRequiredInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = repca.size() + bus.size();
          for (size_t port = 0; port < E::MAXIMUM; ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = repca.getSignals();
          signals.template attachSignalNode<Ext::IBRANCHR>(&input_nodes_[E::IBRANCHR]);
          signals.template attachSignalNode<Ext::IBRANCHI>(&input_nodes_[E::IBRANCHI]);
          signals.template attachSignalNode<Ext::PBRANCH>(&input_nodes_[E::PBRANCH]);
          signals.template attachSignalNode<Ext::QBRANCH>(&input_nodes_[E::QBRANCH]);
          signals.template attachSignalNode<Ext::FREQ>(&input_nodes_[E::FREQ]);
        }

        void attachAllInputs(RealT initial_value = 0.0)
        {
          attachRequiredInputs(initial_value);

          auto& signals = repca.getSignals();
          signals.template attachSignalNode<Ext::FREQREF>(&input_nodes_[E::FREQREF]);
          signals.template attachSignalNode<Ext::VREF>(&input_nodes_[E::VREF]);
          signals.template attachSignalNode<Ext::QREF>(&input_nodes_[E::QREF]);
          signals.template attachSignalNode<Ext::PPLANTREF>(&input_nodes_[E::PPLANTREF]);
        }

        void seedCommands(RealT qext, RealT pext)
        {
          qext_node_.init(static_cast<T>(qext));
          pext_node_.init(static_cast<T>(pext));
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
          seedCommands(qext, pext);
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
          return qext_node_.read();
        }

        T pext() const
        {
          return pext_node_.read();
        }

        T& input(size_t port)
        {
          return input_values_[port];
        }

        IdxT inputIndex(size_t port) const
        {
          return input_indices_[port];
        }

        PhasorDynamics::Bus<T, IdxT>              bus;
        PhasorDynamics::Converter::Repca<T, IdxT> repca;
      };

      static constexpr RealT kStateVr = 0.9;
      static constexpr RealT kStateVi = 0.4;

      static constexpr size_t kBusVrColumn        = I::MAXIMUM;
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

      Data makeDynamicData() const
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
        auto data                         = makeDynamicData();
        data.parameters[Params::dbdupper] = 0.02;
        data.parameters[Params::emin]     = -0.8;
        data.parameters[Params::femin]    = -0.6;
        return data;
      }

      template <typename T>
      void setInitializationInputs(Fixture<T>& fixture) const
      {
        fixture.input(E::IBRANCHR) = static_cast<T>(0.2);
        fixture.input(E::IBRANCHI) = static_cast<T>(-0.1);
        fixture.input(E::PBRANCH)  = static_cast<T>(0.4);
        fixture.input(E::QBRANCH)  = static_cast<T>(0.1);
        fixture.input(E::FREQ)     = static_cast<T>(0.99);
      }

      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(E::IBRANCHR)  = static_cast<T>(0.08);
        fixture.input(E::IBRANCHI)  = static_cast<T>(-0.02);
        fixture.input(E::PBRANCH)   = static_cast<T>(0.41);
        fixture.input(E::QBRANCH)   = static_cast<T>(0.13);
        fixture.input(E::FREQ)      = static_cast<T>(0.99);
        fixture.input(E::FREQREF)   = static_cast<T>(1.0);
        fixture.input(E::VREF)      = static_cast<T>(1.01);
        fixture.input(E::QREF)      = static_cast<T>(0.12);
        fixture.input(E::PPLANTREF) = static_cast<T>(0.55);
      }

      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Converter::Repca<T, IdxT>& repca) const
      {
        setState(repca,
                 {{I::VMEAS, 0.98},
                  {I::QMEAS, 0.11},
                  {I::XQPI, 0.07},
                  {I::XQLAG, 0.14},
                  {I::PMEAS, 0.44},
                  {I::XPPI, 0.21},
                  {I::PREF, 0.60},
                  {I::V, 1.0},
                  {I::VLDC, 0.99},
                  {I::VDROOP, 1.05},
                  {I::VCTRL, 1.02},
                  {I::SFRZ, 0.8},
                  {I::ERQ, 0.03},
                  {I::ERQDB, 0.02},
                  {I::ERQLIM, 0.02},
                  {I::QPI, 0.27},
                  {I::QEXT, 0.20},
                  {I::EF, 0.01},
                  {I::EP, 0.04},
                  {I::EPLIM, 0.04},
                  {I::PPI, 0.66},
                  {I::PEXT, 0.61}});
        setDerivative(repca,
                      {{I::VMEAS, 0.01},
                       {I::QMEAS, -0.02},
                       {I::XQPI, 0.03},
                       {I::XQLAG, -0.04},
                       {I::PMEAS, 0.02},
                       {I::XPPI, -0.01},
                       {I::PREF, 0.05}});
      }

      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), 0.9, 0.4);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), 0.9, 0.4);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        implicit_defaults.input(E::PBRANCH) = 0.2;
        implicit_defaults.input(E::QBRANCH) = 0.1;
        implicit_defaults.input(E::FREQ)    = 1.0;
        explicit_defaults.input(E::PBRANCH) = 0.2;
        explicit_defaults.input(E::QBRANCH) = 0.1;
        explicit_defaults.input(E::FREQ)    = 1.0;

        bool success = implicit_defaults.initialize(0.1, 0.2)
                       && explicit_defaults.initialize(0.1, 0.2);
        if (!success)
        {
          std::cout << "REPCA documented-default comparison failed to initialize\n";
          return false;
        }

        success &= (implicit_defaults.repca.evaluateResidual() == 0);
        success &= (explicit_defaults.repca.evaluateResidual() == 0);
        success &= vectorsMatch(implicit_defaults.repca.y(),
                                explicit_defaults.repca.y(),
                                "documented-default state");
        success &= vectorsMatch(implicit_defaults.repca.yp(),
                                explicit_defaults.repca.yp(),
                                "documented-default derivative");
        success &= vectorsMatch(implicit_defaults.repca.getResidual(),
                                explicit_defaults.repca.getResidual(),
                                "documented-default residual");
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          success &= rowMatches(implicit_defaults.input(port),
                                explicit_defaults.input(port),
                                "documented-default signal",
                                port,
                                "");
        }

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.repca);
        setAnswerKeyState(explicit_defaults.repca);
        success &= (implicit_defaults.repca.evaluateResidual() == 0);
        success &= (explicit_defaults.repca.evaluateResidual() == 0);
        success &= vectorsMatch(implicit_defaults.repca.getResidual(),
                                explicit_defaults.repca.getResidual(),
                                "documented-default dynamic residual");
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

      bool initializationRejectedAtomically(const Data&     data,
                                            RealT           qext,
                                            RealT           pext,
                                            const char*     label,
                                            NonfiniteTarget target = NonfiniteTarget::NONE) const
      {
        Fixture<ScalarT> fixture(data, 0.8, 0.6);
        fixture.attachAllInputs(77.0);
        setInitializationInputs(fixture);
        if (!fixture.prepare(qext, pext))
        {
          return false;
        }

        const RealT infinity = std::numeric_limits<RealT>::infinity();
        if (target == NonfiniteTarget::FREQUENCY)
        {
          fixture.input(E::FREQ) = infinity;
        }
        if (target == NonfiniteTarget::BUS_VOLTAGE)
        {
          fixture.bus.Vr() = infinity;
          fixture.bus.y().setDataUpdated();
        }

        auto* y  = fixture.repca.y().getData();
        auto* yp = fixture.repca.yp().getData();
        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          y[row]  = 0.125 + 0.01 * static_cast<RealT>(row);
          yp[row] = -0.25 - 0.01 * static_cast<RealT>(row);
        }
        // Restore assigned commands after filling model storage with sentinels.
        fixture.seedCommands(qext, pext);
        fixture.repca.y().setDataUpdated();
        fixture.repca.yp().setDataUpdated();

        const auto                    y_before   = copyVector(fixture.repca.y());
        const auto                    yp_before  = copyVector(fixture.repca.yp());
        const auto                    bus_before = copyVector(fixture.bus.y());
        std::array<RealT, E::MAXIMUM> inputs_before{};
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          inputs_before[port] = fixture.input(port);
        }

        bool success = true;
        if (fixture.repca.initialize() == 0)
        {
          std::cout << "Expected REPCA initialization rejection: " << label << '\n';
          success = false;
        }

        success &= valueUnchanged(fixture.qext(), qext, "qext signal", I::QEXT);
        success &= valueUnchanged(fixture.pext(), pext, "pext signal", I::PEXT);
        success &= vectorUnchanged(fixture.repca.y(), y_before, "state");
        success &= vectorUnchanged(fixture.repca.yp(), yp_before, "derivative");
        success &= vectorUnchanged(fixture.bus.y(), bus_before, "bus state");
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          success &= valueUnchanged(fixture.input(port),
                                    inputs_before[port],
                                    "external signal",
                                    port);
        }
        return success;
      }

      template <typename T>
      void setState(PhasorDynamics::Converter::Repca<T, IdxT>& repca, Rows rows) const
      {
        auto* y = repca.y().getData();
        for (const auto& [row, value] : rows)
        {
          y[row] = static_cast<T>(value);
        }
        repca.y().setDataUpdated();
      }

      template <typename T>
      void setDerivative(PhasorDynamics::Converter::Repca<T, IdxT>& repca, Rows rows) const
      {
        auto* yp = repca.yp().getData();
        for (const auto& [row, value] : rows)
        {
          yp[row] = static_cast<T>(value);
        }
        repca.yp().setDataUpdated();
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
        std::cout << " mismatch: " << std::setprecision(16) << actual
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
        std::cout << label << " mismatch: " << std::setprecision(16) << actual
                  << " != " << expected << '\n';
        return false;
      }

      static bool valueUnchanged(RealT       actual,
                                 RealT       expected,
                                 const char* what,
                                 size_t      index)
      {
        if (actual == expected)
        {
          return true;
        }
        std::cout << "REPCA " << what << ' ' << index
                  << " changed: " << std::setprecision(16) << actual
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
          const auto& [row, expected] = rows[i];

          success &= rowMatches(static_cast<RealT>(values[row]),
                                expected,
                                what,
                                row,
                                context);
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

      template <size_t size>
      bool residualsMatch(const RepcaT&                repca,
                          const std::array<Row, size>& rows,
                          const char*                  context = "") const
      {
        return rowsMatch(repca.getResidual(), rows.data(), size, "residual", context);
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
        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          success &= rowMatches(f[row], 0.0, "residual", row, "at rest");
          success &= rowMatches(yp[row], 0.0, "derivative", row, "at rest");
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
          success &= valueUnchanged(static_cast<RealT>(values[row]),
                                    snapshot[row],
                                    what,
                                    row);
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
          success &= rowMatches(static_cast<RealT>(left_values[row]),
                                static_cast<RealT>(right_values[row]),
                                what,
                                row,
                                "");
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

      std::array<JacobianRow, I::MAXIMUM> expectedJacobian() const
      {
        return {{
            {{I::VMEAS, -6.0}, {I::VCTRL, 5.0}},
            {{I::QMEAS, -6.0}, {externalColumn(E::QBRANCH), 10.0}},
            {{I::XQPI, -1.0}, {I::SFRZ, 0.06}, {I::ERQLIM, 2.4}},
            {{I::XQLAG, -1.6666666666666667}, {I::QPI, 0.6666666666666666}},
            {{I::PMEAS, -3.5}, {externalColumn(E::PBRANCH), 5.0}},
            {{I::XPPI, -1.0}, {I::EPLIM, 1.8}},
            {{I::PREF, -3.0}, {I::PPI, 2.0}},
            {{I::V, -2.0}, {kBusVrColumn, 1.8}, {kBusViColumn, 0.8}},
            {{I::VLDC, -1.98},
             {kBusVrColumn, 1.7956},
             {kBusViColumn, 0.796},
             {externalColumn(E::IBRANCHR), -0.059792},
             {externalColumn(E::IBRANCHI), 0.037948}},
            {{I::V, 1.0}, {I::VDROOP, -1.0}, {externalColumn(E::QBRANCH), 0.8}},
            {{I::VLDC, 1.0}, {I::VCTRL, -1.0}},
            {{I::SFRZ, -1.0}},
            {{I::VMEAS, -1.0}, {I::ERQ, -1.0}, {externalColumn(E::VREF), 1.0}},
            {{I::ERQ, 0.50000614417460221}, {I::ERQDB, -1.0}},
            {{I::ERQDB, 1.0}, {I::ERQLIM, -1.0}},
            {{I::XQPI, 1.0}, {I::ERQLIM, 2.0}, {I::QPI, -1.0}},
            {{I::XQLAG, 1.3}, {I::QPI, 0.2}, {I::QEXT, -3.0}},
            {{I::EF, -1.0},
             {externalColumn(E::FREQ), -0.23963778765414229},
             {externalColumn(E::FREQREF), 0.23963778765414229}},
            {{I::PMEAS, -1.0},
             {I::EF, 1.9168273035060777},
             {I::EP, -1.0},
             {externalColumn(E::PPLANTREF), 2.0}},
            {{I::EP, 1.0}, {I::EPLIM, -1.0}},
            {{I::XPPI, 1.0}, {I::EPLIM, 1.7}, {I::PPI, -1.0}},
            {{I::PREF, 1.0}, {I::PEXT, -2.0}},
        }};
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

      bool jacobianMatches(const std::vector<JacobianRow>&            actual,
                           const std::array<JacobianRow, I::MAXIMUM>& expected,
                           const char*                                source,
                           RealT                                      tolerance) const
      {
        if (actual.size() != expected.size())
        {
          std::cout << "REPCA " << source << " Jacobian row-count mismatch\n";
          return false;
        }

        constexpr size_t column_count = externalColumn(E::MAXIMUM);
        bool             success      = true;
        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          for (size_t column = 0; column < column_count; ++column)
          {
            const RealT actual_value   = mapValue(actual[row], column);
            const RealT expected_value = mapValue(expected[row], column);
            if (!isEqual(actual_value, expected_value, tolerance))
            {
              std::cout << "REPCA " << source << " Jacobian row " << row
                        << " column " << column << " mismatch: "
                        << std::setprecision(16) << actual_value << " != "
                        << expected_value << '\n';
              success = false;
            }
          }
          for (const auto& [column, value] : actual[row])
          {
            if (column >= column_count && !isEqual(value, 0.0, tolerance))
            {
              std::cout << "REPCA " << source << " Jacobian row " << row
                        << " has unexpected column " << column << '\n';
              success = false;
            }
          }
        }
        return success;
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y     = fixture.repca.y().getData();
        auto* yp    = fixture.repca.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          y[row].setVariableNumber(row);
          yp[row].setVariableNumber(row);
        }
        for (size_t row = 0; row < static_cast<size_t>(fixture.bus.size()); ++row)
        {
          bus_y[row].setVariableNumber(kBusVrColumn + row);
        }
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          fixture.input(port).setVariableNumber(fixture.inputIndex(port));
        }

        fixture.repca.y().setDataUpdated();
        fixture.repca.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<JacobianRow> dependencyTrackingJacobian(const Data& data,
                                                          TestStatus& success) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.0, 0.0);
        setAnswerKeyState(fixture.repca);
        numberVariables(fixture);
        success *= (fixture.repca.evaluateResidual() == 0);

        std::vector<JacobianRow> rows(I::MAXIMUM);
        const auto*              f = fixture.repca.getResidual().getData();
        for (size_t row = 0; row < I::MAXIMUM; ++row)
        {
          rows[row] = f[row].getDependencies();
        }
        return rows;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      std::vector<JacobianRow> enzymeJacobian(const Data& data,
                                              TestStatus& success) const
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
        fixture.repca.updateTime(0.0, 1.0);
        success *= (fixture.repca.evaluateResidual() == 0);
        success *= (fixture.repca.evaluateJacobian() == 0);
        success *= (fixture.repca.constructCsr() == 0);
        return MapFromCsr(fixture.repca.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
