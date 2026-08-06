#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
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
    class GovernorGastPtiTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      GovernorGastPtiTests()  = default;
      ~GovernorGastPtiTests() = default;

      // Smooth gates, base conversion, and time-constant scaling accumulate a
      // few floating-point operations beyond one machine epsilon.
      static constexpr RealT kTol =
          static_cast<RealT>(10.0) * std::numeric_limits<RealT>::epsilon();

      /// Construction, parameter types and domains, lifecycle, and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing GASTPTI defaults and invalid configurations. "
                         "Logged errors and time-constant warnings are expected.");

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> empty;
        success *= (empty.size() == static_cast<IdxT>(index(Internal::MAXIMUM)));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() > 0); // required pmech assignment is absent

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> configured(makeData());
        success *= (configured.size() == static_cast<IdxT>(index(Internal::MAXIMUM)));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() > 0); // required pmech assignment is absent

        // Framework binding precedes model allocation; verification must not
        // inspect index maps until allocate() has sized them.
        typename GastPtiT::VectorT bound_y;
        typename GastPtiT::VectorT bound_yp;
        typename GastPtiT::VectorT bound_f;
        typename GastPtiT::VectorT bound_abs_tol;
        const auto                 bound_size = static_cast<IdxT>(index(Internal::MAXIMUM));
        bound_y.resize(bound_size);
        bound_yp.resize(bound_size);
        bound_f.resize(bound_size);
        bound_abs_tol.resize(bound_size);
        PhasorDynamics::SignalNode<ScalarT, IdxT> bound_pmech;
        GastPtiT                                  bound(makeData());
        bound.getSignals().template assignSignalNode<Internal::PMECH>(&bound_pmech);
        success *= (bound.bind(bound_y, bound_yp, bound_f, bound_abs_tol, 0) == 0);
        success *= (bound.verify() == 0);

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> minimal(makeMinimalData());
        success *= (minimal.verify() > 0); // required pmech assignment is absent
        success *= (verifyData(makeData()) == 0);
        success *= defaultsMatchDocumentedValues();
        success *= omittedRatingUsesSystemBase();

        success *= invalidParameterCase(Params::R, 0.0);
        success *= invalidParameterCase(Params::R, -0.1);
        success *= invalidParameterCase(Params::T1, -0.1);
        success *= invalidParameterCase(Params::T2, -0.1);
        success *= invalidParameterCase(Params::T3, -0.1);
        success *= invalidParameterCase(Params::At, -0.1);
        success *= invalidParameterCase(Params::Kt, -0.1);
        success *= invalidParameterCase(Params::Vmin, 2.0); // above Vmax
        success *= invalidParameterCase(Params::Dturb, -0.1);
        success *= invalidParameterCase(Params::Trate, 0.0);
        success *= invalidParameterCase(Params::Trate, -1.0);

        const std::array<Params, 10> real_parameters{{
            Params::R,
            Params::T1,
            Params::T2,
            Params::T3,
            Params::At,
            Params::Kt,
            Params::Vmax,
            Params::Vmin,
            Params::Dturb,
            Params::Trate,
        }};
        for (const Params parameter : real_parameters)
        {
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::infinity());
        }

        // A finite serialized rating that overflows its MW-to-VA conversion is
        // still an invalid component base.
        success *= invalidParameterCase(Params::Trate,
                                        std::numeric_limits<RealT>::max());
        success *= invalidParameterCase(Params::Trate,
                                        std::numeric_limits<RealT>::denorm_min());

        // Equal configured limits are valid; reversed limits are not.
        auto equal                      = makeData();
        equal.parameters[Params::Vmin]  = 0.5;
        equal.parameters[Params::Vmax]  = 0.5;
        success                        *= (verifyData(equal) == 0);

        auto reversed                      = makeData();
        reversed.parameters[Params::Vmin]  = 0.6;
        reversed.parameters[Params::Vmax]  = 0.5;
        success                           *= (verifyData(reversed) > 0);

        // Narrow configured limits remain valid.
        auto narrow                      = makeData();
        narrow.parameters[Params::Vmax]  = 0.01;
        success                         *= (verifyData(narrow) == 0);

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                      = makeData();
        integer_real.parameters[Params::Trate] = static_cast<IdxT>(50);
        Fixture<ScalarT> integer_fixture(integer_real);
        success *= integer_fixture.initialize(0.4);
        success *= stateMatches(integer_fixture.gastpti,
                                {{Internal::XFLOW, 0.8}},
                                "integer-valued component base");

        success *= invalidParameterCase(Params::T1, true);

        const std::array<RealT, 5> invalid_system_bases{{
            -100.0e6,
            ZERO<RealT>,
            std::numeric_limits<RealT>::denorm_min(),
            std::numeric_limits<RealT>::quiet_NaN(),
            std::numeric_limits<RealT>::infinity(),
        }};
        for (const RealT system_base : invalid_system_bases)
        {
          Fixture<ScalarT> invalid_base(makeData(), system_base);
          success *= (invalid_base.gastpti.allocate() == 0);
          success *= (invalid_base.gastpti.verify() > 0);
        }

        success *= unlinkedSignalRejected<External::OMEGA>();
        success *= unlinkedSignalRejected<External::PREF>();
        success *= aliasedSignalsRejected();

        Fixture<ScalarT> unallocated(makeData());
        success *= (unallocated.gastpti.initialize() != 0);

        struct TimeConstantCase
        {
          RealT value;
          RealT expected_residual;
        };

        const std::array<TimeConstantCase, 4> time_constant_cases{{
            {0.0, 1.0},
            {0.0005, 1.0},
            {0.001, 1.0},
            {0.002, 0.5},
        }};
        for (const auto& test_case : time_constant_cases)
        {
          auto time_data                   = makeData();
          time_data.parameters[Params::T1] = test_case.value;
          time_data.parameters[Params::T2] = test_case.value;
          time_data.parameters[Params::T3] = test_case.value;

          Fixture<ScalarT> time_fixture(time_data);
          success *= time_fixture.initialize(0.4);
          setState(time_fixture.gastpti,
                   {{Internal::XVALVE, 0.401},
                    {Internal::XFLOW, 0.4},
                    {Internal::XTEMP, 0.399},
                    {Internal::VLV, 0.402}});
          success *= (time_fixture.evaluate() == 0);
          success *= residualsMatch(
              time_fixture.gastpti,
              {{Internal::XVALVE, test_case.expected_residual},
               {Internal::XFLOW, test_case.expected_residual},
               {Internal::XTEMP, test_case.expected_residual}},
              "in-place time-constant floor boundary");
        }

        return success.report(__func__);
      }

      /// Nonidentity base conversion preserves known inputs, publishes the
      /// unknown reference, and initializes signals, monitors, and tags.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        auto data                      = makeData();
        data.parameters[Params::Trate] = 50.0;

        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        fixture.input(index(External::PREF))  = 99.0; // stale value the publication must replace
        success                              *= fixture.initialize(0.4);
        success                              *= (fixture.gastpti.tagDifferentiable() == 0);
        success                              *= (fixture.evaluate() == 0);

        const auto* y  = fixture.gastpti.y().getData();
        success       *= scalarMatches(y[index(Internal::XVALVE)], 0.8, "XVALVE on component base");
        success       *= scalarMatches(y[index(Internal::XFLOW)], 0.8, "XFLOW on component base");
        success       *= scalarMatches(y[index(Internal::XTEMP)], 0.8, "XTEMP on component base");
        success       *= scalarMatches(y[index(Internal::VLOAD)], 0.8, "VLOAD behind the LV gate");
        success       *= scalarMatches(y[index(Internal::VTEMP)], 2.36, "VTEMP at the temperature limit");
        success       *= scalarMatches(y[index(Internal::VLV)], 0.8, "VLV at the fuel flow");
        success       *= scalarPreserved(fixture.pmech(), 0.4, "preserved pmech seed");

        success *= scalarPreserved(fixture.input(index(External::OMEGA)),
                                   0.0,
                                   "preserved omega input");
        success *= scalarMatches(fixture.input(index(External::PREF)), 0.4, "published pref");

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.gastpti.getMonitor());
        std::stringstream monitor_output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, monitor_output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string monitor_header;
        std::string monitor_values;
        std::getline(monitor_output, monitor_header);
        std::getline(monitor_output, monitor_values);
        success              *= (monitor_header == "t,GastPti_gastpti_test_pmech,"
                                                   "GastPti_gastpti_test_xvalve,"
                                                   "GastPti_gastpti_test_xflow,"
                                                   "GastPti_gastpti_test_xtemp,"
                                                   "GastPti_gastpti_test_vload,"
                                                   "GastPti_gastpti_test_vtemp");
        const auto monitored  = Tokenizer<RealT>(monitor_values, ',')();
        if (monitored.size() == 7)
        {
          success *= scalarMatches(monitored[1], 0.4, "monitored pmech");
          success *= scalarMatches(monitored[2], 0.8, "monitored xvalve");
          success *= scalarMatches(monitored[3], 0.8, "monitored xflow");
          success *= scalarMatches(monitored[4], 0.8, "monitored xtemp");
          success *= scalarMatches(monitored[5], 0.8, "monitored vload");
          success *= scalarMatches(monitored[6], 2.36, "monitored vtemp");
        }
        else
        {
          std::cout << "GASTPTI monitor emitted " << monitored.size()
                    << " values instead of 7\n";
          success = false;
        }

        for (size_t i = 0; i < static_cast<size_t>(fixture.gastpti.size()); ++i)
        {
          const bool expected = i <= index(Internal::XTEMP);
          if (fixture.gastpti.tag()[i] != expected)
          {
            std::cout << "GASTPTI differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= allResidualsZero(fixture.gastpti);

        // A system-base reference step lands on the droop row scaled by the
        // base ratio.
        fixture.input(index(External::PREF))  = 0.5; // the published 0.4 plus a 0.1 step
        success                              *= (fixture.evaluate() == 0);
        success                              *= residualsMatch(fixture.gastpti,
                                                               {{Internal::VLOAD, 0.01}},
                                  "reference step on the component base");

        // Unattached ports fall back to the reference latched by
        // initialize(), so the same steady state holds without a controller.
        Fixture<ScalarT> fallback(data);
        success *= fallback.initialize(0.4);
        success *= (fallback.evaluate() == 0);
        success *= allResidualsZero(fallback.gastpti);

        struct SpeedInitializationCase
        {
          RealT omega;
          RealT xflow;
          RealT vtemp;
          RealT pref;
        };

        const std::array<SpeedInitializationCase, 2> speed_cases{{
            {0.05, 0.805, 2.3585, 0.9025},
            {-0.05, 0.795, 2.3615, -0.1025},
        }};
        for (const auto& test_case : speed_cases)
        {
          Fixture<ScalarT> speed_fixture(data);
          speed_fixture.attachAllInputs();
          speed_fixture.input(index(External::OMEGA))  = test_case.omega;
          success                                     *= speed_fixture.initialize(0.4);
          success                                     *= stateMatches(speed_fixture.gastpti,
                                                                      {{Internal::XVALVE, test_case.xflow},
                                                                       {Internal::XFLOW, test_case.xflow},
                                                                       {Internal::XTEMP, test_case.xflow},
                                                                       {Internal::VLOAD, test_case.xflow},
                                                                       {Internal::VTEMP, test_case.vtemp},
                                                                       {Internal::VLV, test_case.xflow}},
                                  "signed nonzero-speed initialization");
          success                                     *= scalarPreserved(speed_fixture.pmech(),
                                     0.4,
                                     "signed-speed pmech preservation");
          success                                     *= scalarPreserved(speed_fixture.input(index(External::OMEGA)),
                                     test_case.omega,
                                     "signed-speed input preservation");
          success                                     *= scalarMatches(speed_fixture.input(index(External::PREF)),
                                   test_case.pref,
                                   "signed-speed pref publication");
          success                                     *= (speed_fixture.evaluate() == 0);
          success                                     *= allResidualsZero(speed_fixture.gastpti);
        }

        return success.report(__func__);
      }

      /// Initialization-domain boundaries, effective limits, and exact failure
      /// atomicity.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible GASTPTI temperature-gate and "
                         "configuration initialization points. Logged errors "
                         "and response-limit warnings are expected.");

        struct RejectionCase
        {
          const char* label;
          RealT       at;
          RealT       pmech;
        };

        // makeResidualData() halves the power base, so a 0.4 seed is a 0.8
        // component-base fuel flow.
        const std::array<RejectionCase, 2> rejected{{
            {"temperature-gate margin at equality", 0.8, 0.4},
            {"temperature-gate margin negative", 0.0, 0.4},
        }};

        for (const auto& test_case : rejected)
        {
          auto data                    = makeResidualData();
          data.parameters[Params::At]  = test_case.at;
          success                     *= initializationRejectedAtomically(
              data, test_case.pmech, test_case.label);
        }

        // An invalid configuration is rejected before any state is written.
        auto invalid_data                  = makeResidualData();
        invalid_data.parameters[Params::R] = 0.0;
        Fixture<ScalarT> invalid_fixture(invalid_data);
        invalid_fixture.attachAllInputs();
        success *= (invalid_fixture.gastpti.allocate() == 0);
        poisonState(invalid_fixture, 0.4);
        const auto invalid_y  = copyVector(invalid_fixture.gastpti.y());
        const auto invalid_yp = copyVector(invalid_fixture.gastpti.yp());
        if (invalid_fixture.gastpti.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: invalid configuration\n";
          success = false;
        }
        success *= vectorUnchanged(invalid_fixture.gastpti.y(), invalid_y, "state");
        success *= vectorUnchanged(invalid_fixture.gastpti.yp(), invalid_yp, "derivative");

        const std::array<RealT, 3> nonfinite_values{{
            std::numeric_limits<RealT>::quiet_NaN(),
            std::numeric_limits<RealT>::infinity(),
            -std::numeric_limits<RealT>::infinity(),
        }};
        for (const RealT pmech : nonfinite_values)
        {
          success *= initializationRejectedAtomically(
              makeResidualData(), pmech, "non-finite pmech seed");
        }

        for (const RealT omega : nonfinite_values)
        {
          success *= initializationRejectedAtomically(
              makeResidualData(), 0.4, "non-finite speed seed", omega);
        }

        // Each case starts from finite parameters and finite seeds, then
        // overflows a different initialization candidate. None may commit.
        auto power_candidate  = makeResidualData();
        success              *= initializationRejectedAtomically(
            power_candidate,
            std::numeric_limits<RealT>::max(),
            "non-finite component-base power candidate");

        auto flow_candidate                       = makeResidualData();
        flow_candidate.parameters[Params::Dturb]  = std::numeric_limits<RealT>::max();
        success                                  *= initializationRejectedAtomically(
            flow_candidate, 0.4, "non-finite fuel-flow candidate", 2.0);

        auto temperature_candidate                    = makeResidualData();
        temperature_candidate.parameters[Params::At]  = std::numeric_limits<RealT>::max();
        success                                      *= initializationRejectedAtomically(
            temperature_candidate, 0.4, "non-finite temperature candidate");

        auto reference_candidate                       = makeResidualData();
        reference_candidate.parameters[Params::R]      = std::numeric_limits<RealT>::denorm_min();
        reference_candidate.parameters[Params::Dturb]  = ZERO<RealT>;
        success                                       *= initializationRejectedAtomically(
            reference_candidate,
            0.4,
            "non-finite reference candidate",
            std::numeric_limits<RealT>::max());

        // Normal response expands its effective bounds to an over-rated initial
        // flow and remains exactly at rest.
        Fixture<ScalarT> over_rated(makeResidualData());
        over_rated.attachAllInputs();
        success *= over_rated.initialize(0.6); // fuel flow 1.2 above Vmax = 1.1
        success *= stateMatches(over_rated.gastpti,
                                {{Internal::XVALVE, 1.2},
                                 {Internal::XFLOW, 1.2},
                                 {Internal::VLV, 1.2}},
                                "over-rated dispatch");
        success *= scalarPreserved(over_rated.pmech(),
                                   0.6,
                                   "preserved over-rated pmech seed");
        success *= (over_rated.evaluate() == 0);
        success *= allResidualsZero(over_rated.gastpti);

        // A failed reinitialization must preserve the last committed effective
        // limits as well as state, derivatives, and pref.
        Fixture<ScalarT> reinitialize(makeResidualData());
        reinitialize.attachAllInputs();
        success *= reinitialize.initialize(0.6);
        reinitialize.seedPmech(1.0);
        const auto y_before    = copyVector(reinitialize.gastpti.y());
        const auto yp_before   = copyVector(reinitialize.gastpti.yp());
        const auto pref_before = reinitialize.input(index(External::PREF));
        if (reinitialize.gastpti.initialize() == 0)
        {
          std::cout << "Expected failed GASTPTI reinitialization\n";
          success = false;
        }
        success *= vectorUnchanged(reinitialize.gastpti.y(), y_before, "reinitialized state");
        success *= vectorUnchanged(reinitialize.gastpti.yp(), yp_before, "reinitialized derivative");
        success *= scalarPreserved(reinitialize.input(index(External::PREF)),
                                   pref_before,
                                   "reinitialized pref");

        reinitialize.seedPmech(0.6);
        setState(reinitialize.gastpti,
                 {{Internal::XVALVE, 1.2},
                  {Internal::VLV, 1.45}});
        setDerivative(reinitialize.gastpti, {{Internal::XVALVE, 0.0}});
        success *= (reinitialize.evaluate() == 0);
        success *= residualsMatch(reinitialize.gastpti,
                                  {{Internal::XVALVE, 0.35714285714285715}},
                                  "failed reinitialization preserves effective limits");

        struct EqualLimitTemperatureCase
        {
          const char* label;
          RealT       at;
          RealT       vtemp;
          RealT       vlv;
        };

        const std::array<EqualLimitTemperatureCase, 2> equal_limit_temperature_cases{{
            {"equal valve limits at the temperature limit",
             0.8,
             0.8,
             0.797111886747667},
            {"equal valve limits above the temperature limit",
             0.0,
             -0.32,
             -0.32},
        }};
        for (const auto& test_case : equal_limit_temperature_cases)
        {
          auto equal_limit_data                     = makeResidualData();
          equal_limit_data.parameters[Params::At]   = test_case.at;
          equal_limit_data.parameters[Params::Vmin] = 0.8;
          equal_limit_data.parameters[Params::Vmax] = 0.8;

          Fixture<ScalarT> equal_limits(equal_limit_data);
          equal_limits.attachAllInputs();
          success *= equal_limits.initialize(0.4);
          success *= stateMatches(equal_limits.gastpti,
                                  {{Internal::XVALVE, 0.8},
                                   {Internal::XFLOW, 0.8},
                                   {Internal::XTEMP, 0.8},
                                   {Internal::VLOAD, 0.8},
                                   {Internal::VTEMP, test_case.vtemp},
                                   {Internal::VLV, test_case.vlv},
                                   {Internal::PMECH, 0.4}},
                                  test_case.label);
          success *= (equal_limits.evaluate() == 0);
          success *= allResidualsZero(equal_limits.gastpti);
        }

        // An unattached reference retains its last successful latch when a
        // later active reinitialization is rejected.
        Fixture<ScalarT> fallback_reinitialize(makeResidualData());
        success *= fallback_reinitialize.initialize(0.4);
        fallback_reinitialize.seedPmech(1.0);
        success *= (fallback_reinitialize.gastpti.initialize() != 0);
        fallback_reinitialize.seedPmech(0.4);
        success *= (fallback_reinitialize.evaluate() == 0);
        success *= allResidualsZero(fallback_reinitialize.gastpti);

        // A zero mechanical-power seed stays admissible.
        Fixture<ScalarT> zero_seed(makeResidualData());
        zero_seed.attachAllInputs();
        success *= zero_seed.initialize(0.0);
        success *= stateMatches(zero_seed.gastpti,
                                {{Internal::XFLOW, 0.0}, {Internal::VTEMP, 2.52}},
                                "zero seed");
        success *= (zero_seed.evaluate() == 0);
        success *= allResidualsZero(zero_seed.gastpti);

        auto negative_data                     = makeResidualData();
        negative_data.parameters[Params::Vmin] = -1.0;
        Fixture<ScalarT> negative_seed(negative_data);
        negative_seed.attachAllInputs();
        success *= negative_seed.initialize(-0.1);
        success *= stateMatches(negative_seed.gastpti,
                                {{Internal::XFLOW, -0.2},
                                 {Internal::VLOAD, -0.2},
                                 {Internal::PMECH, -0.1}},
                                "negative finite dispatch");
        success *= (negative_seed.evaluate() == 0);
        success *= allResidualsZero(negative_seed.gastpti);

        return success.report(__func__);
      }

      /// A fixed near-closed temperature-gate case proves that initialization
      /// uses the inverse smooth ramp and rests all seven residuals exactly
      /// within the documented behavior tolerance.
      TestOutcome initializationExactness()
      {
        TestStatus success = true;

        auto data                   = makeResidualData();
        data.parameters[Params::At] = 0.8 + 1.0e-4 / 1.4;

        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        success *= stateMatches(fixture.gastpti,
                                {{Internal::VLOAD, 0.8155903227031184},
                                 {Internal::VTEMP, 0.8001},
                                 {Internal::VLV, 0.8}},
                                "near-gate initialization");
        success *= scalarMatches(fixture.input(index(External::PREF)),
                                 0.4077951613515592,
                                 "near-gate published pref");
        success *= (fixture.evaluate() == 0);
        success *= allResidualsZero(fixture.gastpti);

        // A very large but finite temperature margin must not erase the
        // ordinary-sized load demand through catastrophic cancellation.
        // With Kt = 0, the exact initialized load demand remains xF0 = 0.8.
        auto large_margin_data                   = makeResidualData();
        large_margin_data.parameters[Params::At] = 1.0e16;
        large_margin_data.parameters[Params::Kt] = 0.0;
        Fixture<ScalarT> large_margin(large_margin_data);
        large_margin.attachAllInputs();
        success *= large_margin.initialize(0.4);
        success *= stateMatches(large_margin.gastpti,
                                {{Internal::VLOAD, 0.8},
                                 {Internal::VLV, 0.8}},
                                "large finite temperature margin");
        success *= (large_margin.evaluate() == 0);
        success *= allResidualsZero(large_margin.gastpti);

        return success.report(__func__);
      }

      /// A fixed numerical answer key for all seven GASTPTI equations. The
      /// expected values are literals, not a second implementation of GASTPTI.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);

        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.gastpti);
        success *= (fixture.evaluate() == 0);

        // Values are pinned after an independent one-time evaluation of the
        // documented equations at setAnswerKeyState()/setAnswerKeyInputs().
        const std::array<VariableValue, index(Internal::MAXIMUM)> expected{{
            {Internal::XVALVE, 0.24614285714285705},
            {Internal::XFLOW, 0.17755555555555544},
            {Internal::XTEMP, -0.001181818181818159},
            {Internal::VLOAD, -0.0326},
            {Internal::VTEMP, 0.978},
            {Internal::VLV, 0.12},
            {Internal::PMECH, -0.1124},
        }};

        success *= (static_cast<size_t>(fixture.gastpti.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.gastpti, expected);

        return success.report(__func__);
      }

      /// Valve anti-windup, speed/damping signs, and adjusted Normal limits.
      TestOutcome governorControl()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing GASTPTI adjusted response limits. "
                         "Logged response-limit warnings are expected.");

        // Both response limits block outward motion and admit restoring motion.
        struct AntiWindupCase
        {
          const char* label;
          RealT       xvalve;
          RealT       vlv;
          RealT       expected;
        };

        const std::array<AntiWindupCase, 4> antiwindup_cases{{
            {"Vmax blocks an outward valve rate", 1.6, 1.85, 0.0},
            {"Vmin blocks an outward valve rate", -0.45, -0.7, 0.0},
            {"Vmax admits a restoring valve rate", 1.6, 1.35, -0.7142857142857143},
            {"Vmin admits a restoring valve rate", -0.45, -0.2, 0.7142857142857143},
        }};
        for (const auto& test_case : antiwindup_cases)
        {
          Fixture<ScalarT> antiwindup(makeResidualData());
          antiwindup.attachAllInputs();
          success *= antiwindup.initialize(0.4);
          setState(antiwindup.gastpti,
                   {{Internal::XVALVE, test_case.xvalve}, {Internal::VLV, test_case.vlv}});
          setDerivative(antiwindup.gastpti, {{Internal::XVALVE, 0.0}});
          success *= (antiwindup.evaluate() == 0);
          success *= residualsMatch(antiwindup.gastpti,
                                    {{Internal::XVALVE, test_case.expected}},
                                    test_case.label);
        }

        // A speed deviation enters the droop and turbine-damping rows.
        Fixture<ScalarT> speed_step(makeResidualData());
        speed_step.attachAllInputs();
        success                                  *= speed_step.initialize(0.4);
        speed_step.input(index(External::OMEGA))  = 0.05;
        success                                  *= (speed_step.evaluate() == 0);
        success                                  *= residualsMatch(speed_step.gastpti,
                                                                   {{Internal::VLOAD, -0.05},
                                                                    {Internal::PMECH, -0.006}},
                                  "speed deviation in the droop and damping rows");

        // Normal response expands both sides of the configured interval to
        // admit the initialized flow. The derived boundary must be used thereafter.
        struct EffectiveBoundaryCase
        {
          RealT pmech;
          RealT command;
          RealT expected;
        };

        const std::array<EffectiveBoundaryCase, 2> effective_boundary_cases{{
            {0.6, 0.25, 0.35714285714285715},
            {0.0, -0.25, -0.35714285714285715},
        }};
        for (const auto& [pmech, command, expected] : effective_boundary_cases)
        {
          Fixture<ScalarT> response(makeResidualData());
          response.attachAllInputs();
          success              *= response.initialize(pmech);
          const RealT boundary  = 2.0 * pmech;
          setState(response.gastpti,
                   {{Internal::XVALVE, boundary},
                    {Internal::VLV, boundary + command}});
          setDerivative(response.gastpti, {{Internal::XVALVE, 0.0}});
          success *= (response.evaluate() == 0);
          success *= residualsMatch(response.gastpti,
                                    {{Internal::XVALVE, expected}},
                                    "adjusted response boundary");
        }

        // The one-half response at an active upper bound has a large state
        // derivative from the smooth gate at an initialization-adjusted bound.
        success *= effectiveLimitBoundaryJacobian(makeResidualData(),
                                                  0.6,
                                                  1.2,
                                                  "adjusted upper bound");

        return success.report(__func__);
      }

      /// The smooth LV gate on both demand sides and at demand equality, plus
      /// the exhaust-temperature feedback row.
      TestOutcome temperatureLimiting()
      {
        TestStatus success = true;

        // The smooth LV gate with the load demand below, above, and equal to
        // the temperature demand.
        struct GateCase
        {
          const char* label;
          RealT       vload;
          RealT       vtemp;
          RealT       expected;
        };

        const std::array<GateCase, 3> gate_cases{{
            {"the load demand wins the LV gate", 0.3, 1.5, 0.3},
            {"the temperature demand wins the LV gate", 1.5, 0.3, 0.3},
            {"equal demands split the smooth LV gate", 0.9, 0.9, 0.897111886747667},
        }};
        for (const auto& test_case : gate_cases)
        {
          Fixture<ScalarT> gate(makeResidualData());
          gate.attachAllInputs();
          success *= gate.initialize(0.4);
          setState(gate.gastpti,
                   {{Internal::VLOAD, test_case.vload},
                    {Internal::VTEMP, test_case.vtemp},
                    {Internal::VLV, 0.0}});
          success *= (gate.evaluate() == 0);
          success *= residualsMatch(gate.gastpti,
                                    {{Internal::VLV, test_case.expected}},
                                    test_case.label);
        }

        // The exhaust-temperature feedback drives the temperature demand.
        Fixture<ScalarT> feedback(makeResidualData());
        feedback.attachAllInputs();
        success *= feedback.initialize(0.4);
        setState(feedback.gastpti,
                 {{Internal::XTEMP, 0.9}, {Internal::VTEMP, 1.1}});
        success *= (feedback.evaluate() == 0);
        success *= residualsMatch(feedback.gastpti,
                                  {{Internal::VTEMP, 1.06}},
                                  "temperature feedback");

        // At equality, the smooth low-value selector splits its sensitivity
        // evenly between the two demand signals.
        Fixture<DependencyTracking::Variable> selector(makeResidualData());
        selector.attachAllInputs();
        success *= selector.initialize(0.4);
        setState(selector.gastpti,
                 {{Internal::VLOAD, 0.9},
                  {Internal::VTEMP, 0.9},
                  {Internal::VLV, 0.7}});
        numberVariables(selector);
        success *= (selector.evaluate() == 0);

        const DependencyTracking::Variable::DependencyMap expected{
            {index(Internal::VLOAD), 0.5},
            {index(Internal::VTEMP), 0.5},
            {index(Internal::VLV), -1.0},
        };
        success *= jacobianRowMatches(
            selector.gastpti.getResidual().getData()[index(Internal::VLV)].getDependencies(),
            expected,
            index(Internal::VLV),
            "selector equality");

        return success.report(__func__);
      }

      /// Check every DependencyTracking row against independent numerical and
      /// structural answer keys in active and collapsed valve configurations.
      TestOutcome dependencyTracking()
      {
        TestStatus success = true;

        const auto                                                     data     = makeResidualData();
        const auto                                                     interior = dependencyTrackingJacobian(data, success);
        const std::vector<DependencyTracking::Variable::DependencyMap> interior_expected{
            {{index(Internal::XVALVE), -3.8571428571428572},
             {index(Internal::VLV), 2.8571428571428572}},
            {{index(Internal::XVALVE), 2.2222222222222223},
             {index(Internal::XFLOW), -3.2222222222222223}},
            {{index(Internal::XFLOW), 0.45454545454545453},
             {index(Internal::XTEMP), -1.4545454545454546}},
            {{index(Internal::VLOAD), -0.06},
             {index(Internal::MAXIMUM) + index(External::OMEGA), -1.0},
             {index(Internal::MAXIMUM) + index(External::PREF), 0.12}},
            {{index(Internal::XTEMP), -0.4},
             {index(Internal::VTEMP), -1.0}},
            {{index(Internal::VLOAD), 1.0},
             {index(Internal::VTEMP), 0.0},
             {index(Internal::VLV), -1.0}},
            {{index(Internal::XFLOW), 1.0},
             {index(Internal::PMECH), -2.0},
             {index(Internal::MAXIMUM) + index(External::OMEGA), -0.12}},
        };
        success *= jacobianMatches(interior,
                                   interior_expected,
                                   "interior answer key");

        // Equal configured limits collapse the effective interval and retain
        // the explicitly stored zero valve-drive dependency.
        auto equal_limit_data                     = data;
        equal_limit_data.parameters[Params::Vmin] = 0.8;
        equal_limit_data.parameters[Params::Vmax] = 0.8;
        const auto equal_limits                   = dependencyTrackingJacobian(equal_limit_data, success);
        const std::vector<DependencyTracking::Variable::DependencyMap>
            equal_limit_expected{
                {{index(Internal::XVALVE), -1.0},
                 {index(Internal::VLV), 0.0}},
                interior_expected[index(Internal::XFLOW)],
                interior_expected[index(Internal::XTEMP)],
                interior_expected[index(Internal::VLOAD)],
                interior_expected[index(Internal::VTEMP)],
                interior_expected[index(Internal::VLV)],
                interior_expected[index(Internal::PMECH)],
            };
        success *= jacobianMatches(equal_limits,
                                   equal_limit_expected,
                                   "equal configured limits");

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Enzyme rows must match dependency tracking in active and collapsed
      /// configurations, and the raw COO pattern must remain fixed as smooth
      /// gates and initialization masks change on the same model instance.
      TestOutcome jacobian()
      {
        noteExpectedLogs("Testing GASTPTI collapsed-to-active Jacobian structure. "
                         "A logged response-limit warning is expected.");

        TestStatus success = true;

        const auto data = makeResidualData();

        const auto dependency_jacobian = dependencyTrackingJacobian(data, success);
        const auto enzyme_jacobian     = enzymeJacobian(data, success);

        success *= jacobianMatches(enzyme_jacobian,
                                   dependency_jacobian,
                                   "Enzyme versus dependency tracking",
                                   kTol);

        auto collapsed_data                     = data;
        collapsed_data.parameters[Params::Vmin] = 0.8;
        collapsed_data.parameters[Params::Vmax] = 0.8;
        const auto collapsed_dependency =
            dependencyTrackingJacobian(collapsed_data, success);
        const auto collapsed_enzyme  = enzymeJacobian(collapsed_data, success);
        success                     *= jacobianMatches(collapsed_enzyme,
                                   collapsed_dependency,
                                   "collapsed Enzyme versus dependency tracking",
                                   kTol);

        success *= selectorAndAntiwindupCooPatternIsStable(data);
        success *= reinitializedCooPatternIsStable(data);

        return success.report(__func__);
      }
#endif

    private:
      using GastPtiT = PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;
      using Data     = typename GastPtiT::ModelDataT;
      using Params   = typename Data::Parameters;
      using Mon      = typename Data::MonitorableVariables;
      using Internal = typename GastPtiT::InternalVariablesT;
      using External = typename GastPtiT::ExternalVariablesT;

      static constexpr size_t index(Internal variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(External variable)
      {
        return static_cast<size_t>(variable);
      }

      struct VariableValue
      {
        Internal variable;
        RealT    value;
      };

      /// Owns the GASTPTI model, the assigned mechanical-power node, and the
      /// attached input nodes. Signal storage is declared before the model so
      /// every referenced node outlives GASTPTI. Copying would invalidate the
      /// model and signal-node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, index(External::MAXIMUM)>                                   input_values_{};
        std::array<IdxT, index(External::MAXIMUM)>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, index(External::MAXIMUM)> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> pmech_node_;

      public:
        explicit Fixture(const Data& data, RealT system_va_base = 100.0e6)
          : gastpti(data)
        {
          gastpti.setSystemBase(60.0, system_va_base);
          gastpti.getSignals().template assignSignalNode<Internal::PMECH>(&pmech_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Attach fixture-owned storage to every external input.
        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = gastpti.size();

          for (size_t port = 0; port < index(External::MAXIMUM); ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = gastpti.getSignals();
          signals.template attachSignalNode<External::OMEGA>(&input_nodes_[index(External::OMEGA)]);
          signals.template attachSignalNode<External::PREF>(&input_nodes_[index(External::PREF)]);
        }

        /// Seed the assigned mechanical-power node on the system base.
        void seedPmech(RealT pmech)
        {
          pmech_node_.init(static_cast<T>(pmech));
        }

        /// Everything GASTPTI initialization requires: allocation,
        /// verification, and a machine-seeded mechanical-power node.
        bool prepare(RealT pmech)
        {
          const bool success = (gastpti.allocate() == 0) && (gastpti.verify() == 0);
          if (!success)
          {
            std::cout << "GASTPTI fixture preparation failed\n";
            return false;
          }

          seedPmech(pmech);
          return true;
        }

        /// prepare() plus successful GASTPTI initialization.
        bool initialize(RealT pmech)
        {
          if (!prepare(pmech))
          {
            return false;
          }
          if (gastpti.initialize() != 0)
          {
            std::cout << "GASTPTI initialization failed\n";
            return false;
          }
          return true;
        }

        int evaluate()
        {
          return gastpti.evaluateResidual();
        }

        T pmech() const
        {
          return pmech_node_.read();
        }

        T& input(size_t port)
        {
          return input_values_[port];
        }

        IdxT inputIndex(size_t port) const
        {
          return input_indices_[port];
        }

        PhasorDynamics::Governor::GastPti<T, IdxT> gastpti;
      };

      Data makeMinimalData() const
      {
        Data data;
        data.device_class          = "GastPti";
        data.disambiguation_string = "gastpti_test";
        data.monitored_variables.insert(Mon::pmech);
        data.monitored_variables.insert(Mon::xvalve);
        data.monitored_variables.insert(Mon::xflow);
        data.monitored_variables.insert(Mon::xtemp);
        data.monitored_variables.insert(Mon::vload);
        data.monitored_variables.insert(Mon::vtemp);
        return data;
      }

      Data makeExplicitDefaultData() const
      {
        auto data = makeMinimalData();

        // These are the documented defaults, spelled out parameter by
        // parameter.
        data.parameters[Params::R]     = 0.05;
        data.parameters[Params::T1]    = 0.4;
        data.parameters[Params::T2]    = 0.1;
        data.parameters[Params::T3]    = 3.0;
        data.parameters[Params::At]    = 1.0;
        data.parameters[Params::Kt]    = 2.0;
        data.parameters[Params::Vmax]  = 1.0;
        data.parameters[Params::Vmin]  = 0.0;
        data.parameters[Params::Dturb] = 0.0;
        data.parameters[Params::Trate] = 100.0;
        return data;
      }

      Data makeData() const
      {
        auto data = makeMinimalData();

        // Finite nondefault values used by routine fixtures.
        data.parameters[Params::R]     = 0.05;
        data.parameters[Params::T1]    = 0.4;
        data.parameters[Params::T2]    = 0.5;
        data.parameters[Params::T3]    = 0.25;
        data.parameters[Params::At]    = 2.0;
        data.parameters[Params::Kt]    = 0.3;
        data.parameters[Params::Vmax]  = 1.2;
        data.parameters[Params::Vmin]  = 0.0;
        data.parameters[Params::Dturb] = 0.1;
        data.parameters[Params::Trate] = 100.0;
        return data;
      }

      Data makeResidualData() const
      {
        auto data = makeData();

        // Dynamic-response parameters: every gain, lag, and limit is
        // nontrivial and the power bases differ.
        data.parameters[Params::Trate] = 50.0;
        data.parameters[Params::R]     = 0.06;
        data.parameters[Params::T1]    = 0.35;
        data.parameters[Params::T2]    = 0.45;
        data.parameters[Params::T3]    = 2.2;
        data.parameters[Params::At]    = 1.8;
        data.parameters[Params::Kt]    = 0.4;
        data.parameters[Params::Vmax]  = 1.1;
        data.parameters[Params::Vmin]  = 0.05;
        data.parameters[Params::Dturb] = 0.12;
        return data;
      }

      /// The external inputs the residual answer key is evaluated against.
      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(index(External::OMEGA)) = static_cast<T>(0.02);
        fixture.input(index(External::PREF))  = static_cast<T>(0.31);
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti) const
      {
        setState(gastpti,
                 {{Internal::XVALVE, 0.62},
                  {Internal::XFLOW, 0.55},
                  {Internal::XTEMP, 0.48},
                  {Internal::VLOAD, 0.83},
                  {Internal::VTEMP, 1.35},
                  {Internal::VLV, 0.71},
                  {Internal::PMECH, 0.33}});
        setDerivative(gastpti,
                      {{Internal::XVALVE, 0.011},
                       {Internal::XFLOW, -0.022},
                       {Internal::XTEMP, 0.033}});
      }

      /// Omitting every parameter must give exactly the model built from the
      /// defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData());
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData());
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        bool success = implicit_defaults.initialize(0.3)
                       && explicit_defaults.initialize(0.3);
        if (!success)
        {
          std::cout << "GASTPTI documented-default comparison failed to initialize\n";
          return false;
        }

        if (implicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (explicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.gastpti.y(),
                             copyVector(explicit_defaults.gastpti.y()),
                             "documented-default state"))
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.gastpti.yp(),
                             copyVector(explicit_defaults.gastpti.yp()),
                             "documented-default derivative"))
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.gastpti.getResidual(),
                             copyVector(explicit_defaults.gastpti.getResidual()),
                             "documented-default residual"))
        {
          success = false;
        }

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.gastpti);
        setAnswerKeyState(explicit_defaults.gastpti);
        if (implicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (explicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.gastpti.getResidual(),
                             copyVector(explicit_defaults.gastpti.getResidual()),
                             "documented-default dynamic residual"))
        {
          success = false;
        }
        return success;
      }

      /// An omitted rating follows a non-100-MVA system base exactly.
      bool omittedRatingUsesSystemBase() const
      {
        constexpr RealT system_va_base = static_cast<RealT>(75.0e6);

        auto omitted_data = makeExplicitDefaultData();
        omitted_data.parameters.erase(Params::Trate);
        auto explicit_data                      = omitted_data;
        explicit_data.parameters[Params::Trate] = static_cast<RealT>(75.0);

        Fixture<ScalarT> omitted(omitted_data, system_va_base);
        Fixture<ScalarT> explicit_system_base(explicit_data, system_va_base);
        omitted.attachAllInputs();
        explicit_system_base.attachAllInputs();

        bool success = omitted.initialize(0.3)
                       && explicit_system_base.initialize(0.3);
        if (!success)
        {
          std::cout << "GASTPTI omitted-rating comparison failed to initialize\n";
          return false;
        }

        if (omitted.evaluate() != 0 || explicit_system_base.evaluate() != 0)
        {
          success = false;
        }
        if (!vectorUnchanged(omitted.gastpti.y(),
                             copyVector(explicit_system_base.gastpti.y()),
                             "omitted-rating state"))
        {
          success = false;
        }
        if (!vectorUnchanged(omitted.gastpti.getResidual(),
                             copyVector(explicit_system_base.gastpti.getResidual()),
                             "omitted-rating residual"))
        {
          success = false;
        }
        if (!scalarMatches(omitted.input(index(External::PREF)),
                           explicit_system_base.input(index(External::PREF)),
                           "omitted-rating published pref"))
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
        return verifyData(data) > 0;
      }

      int verifyData(const Data& data) const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>        pmech;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> model(data);
        model.getSignals().template assignSignalNode<Internal::PMECH>(&pmech);
        return model.verify();
      }

      template <External variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>        unlinked_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>        pmech_node;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> model(makeData());
        model.getSignals().template assignSignalNode<Internal::PMECH>(&pmech_node);
        model.getSignals().template attachSignalNode<variable>(&unlinked_node);
        return model.verify() > 0;
      }

      bool aliasedSignalsRejected() const
      {
        using ModelT = PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;
        using NodeT  = PhasorDynamics::SignalNode<ScalarT, IdxT>;

        bool success = true;

        NodeT  pmech_pref;
        ModelT pref_alias(makeData());
        pref_alias.getSignals().template assignSignalNode<Internal::PMECH>(&pmech_pref);
        pref_alias.getSignals().template attachSignalNode<External::PREF>(&pmech_pref);
        if (pref_alias.allocate() != 0)
        {
          success = false;
        }
        if (!(pref_alias.verify() > 0))
        {
          success = false;
        }

        NodeT  pmech_speed;
        ModelT speed_alias(makeData());
        speed_alias.getSignals().template assignSignalNode<Internal::PMECH>(&pmech_speed);
        speed_alias.getSignals().template attachSignalNode<External::OMEGA>(&pmech_speed);
        if (speed_alias.allocate() != 0)
        {
          success = false;
        }
        if (!(speed_alias.verify() > 0))
        {
          success = false;
        }

        ScalarT shared_value{ZERO<RealT>};
        IdxT    shared_index{static_cast<IdxT>(99)};
        NodeT   shared_input;
        NodeT   pmech;
        shared_input.set(&shared_value, &shared_index);

        ModelT input_alias(makeData());
        input_alias.getSignals().template assignSignalNode<Internal::PMECH>(&pmech);
        input_alias.getSignals().template attachSignalNode<External::OMEGA>(&shared_input);
        input_alias.getSignals().template attachSignalNode<External::PREF>(&shared_input);
        if (input_alias.allocate() != 0)
        {
          success = false;
        }
        if (!(input_alias.verify() > 0))
        {
          success = false;
        }

        return success;
      }

      template <typename VectorT>
      std::vector<RealT> copyVector(const VectorT& vector) const
      {
        const auto* values = vector.getData();
        return std::vector<RealT>(values,
                                  values + static_cast<size_t>(vector.getSize()));
      }

      /// Every row of a vector still holds its snapshot value.
      template <typename VectorT>
      bool vectorUnchanged(const VectorT&            vector,
                           const std::vector<RealT>& snapshot,
                           const char*               what) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t i = 0; i < snapshot.size(); ++i)
        {
          // A rejected initialization may leave a non-finite seed in place.
          const RealT value = static_cast<RealT>(values[i]);
          if (!preserved(value, snapshot[i])
              && !rowMatches(value, snapshot[i], what, i, "changed"))
          {
            success = false;
          }
        }
        return success;
      }

      /// Fill the state and derivative with a recognizable ramp, then re-seed
      /// the aliased pmech entry, so any write by a rejected initialization
      /// is visible.
      void poisonState(Fixture<ScalarT>& fixture, RealT pmech) const
      {
        auto* y  = fixture.gastpti.y().getData();
        auto* yp = fixture.gastpti.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.gastpti.y().getSize()); ++i)
        {
          y[i]  = 0.125 + 0.01 * static_cast<RealT>(i);
          yp[i] = -0.25 - 0.01 * static_cast<RealT>(i);
        }
        fixture.seedPmech(pmech);
        fixture.gastpti.y().setDataUpdated();
        fixture.gastpti.yp().setDataUpdated();
      }

      bool initializationRejectedAtomically(const Data& data,
                                            RealT       pmech,
                                            const char* label,
                                            RealT       omega = ZERO<RealT>) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        fixture.input(index(External::OMEGA)) = omega;
        fixture.input(index(External::PREF))  = 77.0; // must stay untouched on rejection
        if (!fixture.prepare(pmech))
        {
          return false;
        }

        poisonState(fixture, pmech);
        const auto y_before  = copyVector(fixture.gastpti.y());
        const auto yp_before = copyVector(fixture.gastpti.yp());

        bool success = true;
        if (fixture.gastpti.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: " << label << "\n";
          success = false;
        }

        if (!scalarPreserved(fixture.pmech(),
                             pmech,
                             "rejected pmech seed preservation"))
        {
          success = false;
        }
        if (!scalarPreserved(fixture.input(index(External::OMEGA)),
                             omega,
                             "rejected omega preservation"))
        {
          success = false;
        }
        if (!scalarPreserved(fixture.input(index(External::PREF)),
                             77.0,
                             "rejected pref preservation"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.gastpti.y(), y_before, "state"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.gastpti.yp(), yp_before, "derivative"))
        {
          success = false;
        }
        return success;
      }

      /// Write state rows and publish the update, folding in the
      /// setDataUpdated() that a hand-written write block has to remember.
      template <typename T>
      void setState(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti,
                    std::initializer_list<VariableValue>        values) const
      {
        auto* y = gastpti.y().getData();
        for (const auto& [variable, value] : values)
        {
          y[index(variable)] = static_cast<T>(value);
        }
        gastpti.y().setDataUpdated();
      }

      /// setState() for the derivative vector.
      template <typename T>
      void setDerivative(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti,
                         std::initializer_list<VariableValue>        values) const
      {
        auto* yp = gastpti.yp().getData();
        for (const auto& [variable, value] : values)
        {
          yp[index(variable)] = static_cast<T>(value);
        }
        gastpti.yp().setDataUpdated();
      }

      /// Compare one vector row against its expected value. Every row check
      /// in this suite reports through here, so failures share one format.
      /// Variables are named by the canonical enum used by the expectation.
      static bool rowMatches(RealT       actual,
                             RealT       expected,
                             const char* what,
                             size_t      row,
                             const char* context)
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "GASTPTI " << what << " row " << row << ' ' << context
                  << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual
                  << " != " << expected << '\n';
        return false;
      }

      /// Check selected rows of a model vector against expected values.
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
          if (!rowMatches(static_cast<RealT>(vector_values[row]), expected, what, row, context))
          {
            success = false;
          }
        }
        return success;
      }

      bool residualsMatch(const GastPtiT&                      gastpti,
                          std::initializer_list<VariableValue> values,
                          const char*                          context = "") const
      {
        return rowsMatch(gastpti.getResidual(), values, "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const GastPtiT&                        gastpti,
                          const std::array<VariableValue, size>& values,
                          const char*                            context = "") const
      {
        return rowsMatch(gastpti.getResidual(), values, "residual", context);
      }

      bool stateMatches(const GastPtiT&                      gastpti,
                        std::initializer_list<VariableValue> values,
                        const char*                          context = "") const
      {
        return rowsMatch(gastpti.y(), values, "state", context);
      }

      template <size_t size>
      bool stateMatches(const GastPtiT&                        gastpti,
                        const std::array<VariableValue, size>& values,
                        const char*                            context = "") const
      {
        return rowsMatch(gastpti.y(), values, "state", context);
      }

      /// The model sits at a steady state: every residual and every
      /// derivative is zero.
      bool allResidualsZero(const GastPtiT& gastpti) const
      {
        bool        success = true;
        const auto* f       = gastpti.getResidual().getData();
        const auto* yp      = gastpti.yp().getData();
        for (size_t row = 0; row < static_cast<size_t>(gastpti.getResidual().getSize()); ++row)
        {
          if (!rowMatches(static_cast<RealT>(f[row]), 0.0, "residual", row, "at rest"))
          {
            success = false;
          }
          if (!rowMatches(static_cast<RealT>(yp[row]), 0.0, "derivative", row, "at rest"))
          {
            success = false;
          }
        }
        return success;
      }

      bool scalarMatches(ScalarT     actual,
                         ScalarT     expected,
                         const char* label,
                         ScalarT     tolerance = kTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10) << actual
                  << " != " << expected << "\n";
        return false;
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

      bool scalarPreserved(ScalarT actual, ScalarT expected, const char* label) const
      {
        const RealT actual_value   = static_cast<RealT>(actual);
        const RealT expected_value = static_cast<RealT>(expected);
        if (preserved(actual_value, expected_value))
        {
          return true;
        }
        std::cout << label << " changed: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual_value << " != " << expected_value << '\n';
        return false;
      }

      void noteExpectedLogs(const char* message) const
      {
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << "\n";
        Log::setVerbosity(previous_verbosity);
      }

      bool effectiveLimitBoundaryJacobian(const Data& data,
                                          RealT       pmech,
                                          RealT       boundary,
                                          const char* label) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data);
        fixture.attachAllInputs();

        bool success = fixture.initialize(pmech);
        setState(fixture.gastpti,
                 {{Internal::XVALVE, boundary},
                  {Internal::VLV, boundary + 0.25}});
        setDerivative(fixture.gastpti, {{Internal::XVALVE, 0.0}});
        numberVariables(fixture);
        if (fixture.evaluate() != 0)
        {
          success = false;
        }

        const DependencyTracking::Variable::DependencyMap expected{
            {index(Internal::XVALVE), -45.285714285714285},
            {index(Internal::VLV), 1.4285714285714286},
        };
        if (!jacobianRowMatches(
                fixture.gastpti.getResidual().getData()[index(Internal::XVALVE)].getDependencies(),
                expected,
                index(Internal::XVALVE),
                label))
        {
          success = false;
        }
        return success;
      }

      bool jacobianRowMatches(
          const DependencyTracking::Variable::DependencyMap& actual,
          const DependencyTracking::Variable::DependencyMap& expected,
          size_t                                             row,
          const char*                                        context,
          RealT                                              tolerance = kTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }

        std::cout << "GASTPTI Jacobian row " << row << ' ' << context
                  << " mismatch\n";
        return false;
      }

      bool jacobianMatches(
          const std::vector<DependencyTracking::Variable::DependencyMap>& actual,
          const std::vector<DependencyTracking::Variable::DependencyMap>& expected,
          const char*                                                     context,
          RealT                                                           tolerance = kTol) const
      {
        if (actual.size() != expected.size())
        {
          std::cout << "GASTPTI Jacobian " << context << " row-count mismatch: "
                    << actual.size() << " != " << expected.size() << '\n';
          return false;
        }

        bool success = true;
        for (size_t row = 0; row < actual.size(); ++row)
        {
          if (!jacobianRowMatches(actual[row], expected[row], row, context, tolerance))
          {
            success = false;
          }
        }
        return success;
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y  = fixture.gastpti.y().getData();
        auto* yp = fixture.gastpti.yp().getData();

        const auto model_size = static_cast<size_t>(fixture.gastpti.size());
        for (size_t i = 0; i < model_size; ++i)
        {
          y[i].setVariableNumber(i);
          yp[i].setVariableNumber(i);
        }
        for (size_t port = 0; port < index(External::MAXIMUM); ++port)
        {
          fixture.input(port).setVariableNumber(fixture.inputIndex(port));
        }

        fixture.gastpti.y().setDataUpdated();
        fixture.gastpti.yp().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          const Data& data,
          TestStatus& success) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.gastpti);
        numberVariables(fixture);
        success *= (fixture.evaluate() == 0);

        const auto                                               model_size = static_cast<size_t>(fixture.gastpti.size());
        std::vector<DependencyTracking::Variable::DependencyMap> rows(model_size);
        const auto*                                              f = fixture.gastpti.getResidual().getData();
        for (size_t i = 0; i < model_size; ++i)
        {
          rows[i] = f[i].getDependencies();
        }
        return rows;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      using RawCooPattern = std::vector<std::array<IdxT, 2>>;

      RawCooPattern captureRawCooPattern(Fixture<ScalarT>& fixture,
                                         const char*       context) const
      {
        if (fixture.evaluate() != 0)
        {
          std::cout << "GASTPTI residual evaluation failed for " << context << '\n';
          return {};
        }
        if (fixture.gastpti.evaluateJacobian() != 0)
        {
          std::cout << "GASTPTI Jacobian evaluation failed for " << context << '\n';
          return {};
        }

        auto* jacobian = fixture.gastpti.getCooJacobian();
        if (jacobian == nullptr)
        {
          std::cout << "GASTPTI raw COO is absent for " << context << '\n';
          return {};
        }

        constexpr IdxT expected_raw_nnz = 20;
        const IdxT     model_nnz        = fixture.gastpti.nnz();
        const IdxT     coo_nnz          = jacobian->getNnz();
        if (model_nnz != coo_nnz || coo_nnz != expected_raw_nnz)
        {
          std::cout << "GASTPTI raw COO count mismatch for " << context << ": "
                    << model_nnz << " model, " << coo_nnz << " cached, "
                    << expected_raw_nnz << " expected\n";
          return {};
        }

        RawCooPattern pattern;
        pattern.reserve(static_cast<size_t>(coo_nnz));
        const auto* rows    = jacobian->getRowData();
        const auto* columns = jacobian->getColData();
        for (IdxT entry = 0; entry < coo_nnz; ++entry)
        {
          pattern.push_back({rows[entry], columns[entry]});
        }
        return pattern;
      }

      bool selectorAndAntiwindupCooPatternIsStable(const Data& data) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        if (!fixture.initialize(0.4))
        {
          return false;
        }
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.gastpti);
        fixture.gastpti.updateTime(0.0, 1.0);

        setState(fixture.gastpti,
                 {{Internal::VLOAD, 0.9}, {Internal::VTEMP, 0.9}});
        const auto expected = captureRawCooPattern(fixture, "equal selector inputs");
        if (expected.empty())
        {
          return false;
        }

        struct PatternCase
        {
          const char* label;
          Internal    first_variable;
          RealT       first_value;
          Internal    second_variable;
          RealT       second_value;
        };

        const std::array<PatternCase, 4> cases{{
            {"load-limited selector", Internal::VLOAD, 0.3, Internal::VTEMP, 1.5},
            {"temperature-limited selector", Internal::VLOAD, 1.5, Internal::VTEMP, 0.3},
            {"blocked upper valve response", Internal::XVALVE, 1.6, Internal::VLV, 1.85},
            {"restoring upper valve response", Internal::XVALVE, 1.6, Internal::VLV, 1.35},
        }};

        bool success = true;
        for (const auto& test_case : cases)
        {
          setState(fixture.gastpti,
                   {{test_case.first_variable, test_case.first_value},
                    {test_case.second_variable, test_case.second_value}});
          const auto actual = captureRawCooPattern(fixture, test_case.label);
          if (actual != expected)
          {
            std::cout << "GASTPTI raw COO pattern mismatch for "
                      << test_case.label << '\n';
            success = false;
          }
        }
        return success;
      }

      bool reinitializedCooPatternIsStable(const Data& source_data) const
      {
        auto data                     = source_data;
        data.parameters[Params::Vmin] = 0.8;
        data.parameters[Params::Vmax] = 0.8;

        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        if (!fixture.initialize(0.4))
        {
          return false;
        }
        fixture.gastpti.updateTime(0.0, 1.0);

        const auto expected = captureRawCooPattern(fixture, "collapsed valve interval");
        if (expected.empty())
        {
          return false;
        }

        fixture.seedPmech(0.6);
        if (fixture.gastpti.initialize() != 0)
        {
          std::cout << "GASTPTI active-interval reinitialization failed\n";
          return false;
        }

        const auto actual = captureRawCooPattern(
            fixture,
            "active valve interval after reinitialization");
        if (actual != expected)
        {
          std::cout << "GASTPTI raw COO pattern changed across reinitialization\n";
          return false;
        }
        return true;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian(
          const Data& data,
          TestStatus& success) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.gastpti);
        fixture.gastpti.updateTime(0.0, 1.0);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.gastpti.evaluateJacobian() == 0);
        success *= (fixture.gastpti.constructCsr() == 0);
        return MapFromCsr(fixture.gastpti.getCsrJacobian());
      }

#endif
    };
  } // namespace Testing
} // namespace GridKit
