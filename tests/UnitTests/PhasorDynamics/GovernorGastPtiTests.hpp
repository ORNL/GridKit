#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <tuple>
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

      // Behavioral answer keys differ only by accumulated floating-point
      // roundoff.
      static constexpr RealT kBehaviorTol =
          100 * std::numeric_limits<RealT>::epsilon();

      // Enzyme and dependency tracking traverse the same smooth expressions
      // differently; their double-precision derivatives agree to O(1e-10).
      static constexpr RealT kJacobianTol = 1.0e-9;

      /// Construction, parameter types and domains, response modes, lifecycle,
      /// and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        static_assert(static_cast<size_t>(Mode::Normal) == 0);
        static_assert(static_cast<size_t>(Mode::DownOnly) == 1);
        static_assert(static_cast<size_t>(Mode::Fixed) == 2);

        noteExpectedLogs("Testing GASTPTI defaults and invalid configurations. "
                         "Logged errors and time-constant warnings are expected.");

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> empty;
        success *= (empty.size() == static_cast<IdxT>(index(Vars::MAXIMUM)));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() > 0); // required pmech assignment is absent

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> configured(makeData());
        success *= (configured.size() == static_cast<IdxT>(index(Vars::MAXIMUM)));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() > 0); // required pmech assignment is absent

        // Framework binding precedes model allocation; verification must not
        // inspect index maps until allocate() has sized them.
        typename GastPtiT::VectorT bound_y;
        typename GastPtiT::VectorT bound_yp;
        typename GastPtiT::VectorT bound_f;
        typename GastPtiT::VectorT bound_abs_tol;
        const auto                 bound_size = static_cast<IdxT>(index(Vars::MAXIMUM));
        bound_y.resize(bound_size);
        bound_yp.resize(bound_size);
        bound_f.resize(bound_size);
        bound_abs_tol.resize(bound_size);
        PhasorDynamics::SignalNode<ScalarT, IdxT> bound_pmech;
        GastPtiT                                  bound(makeData());
        bound.getSignals().template assignSignalNode<Vars::PMECH>(&bound_pmech);
        success *= (bound.bind(bound_y, bound_yp, bound_f, bound_abs_tol, 0) == 0);
        success *= (bound.verify() == 0);

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> minimal(makeMinimalData());
        success *= (minimal.verify() > 0); // required pmech assignment is absent
        success *= (verifyData(makeData()) == 0);
        success *= defaultsMatchDocumentedValues();

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
          success *= invalidParameterCase(
              parameter, std::numeric_limits<RealT>::quiet_NaN());
        }

        // A finite serialized rating that overflows its MW-to-VA conversion is
        // still an invalid component base.
        success *= invalidParameterCase(Params::Trate,
                                        std::numeric_limits<RealT>::max());
        success *= invalidParameterCase(Params::Trate,
                                        std::numeric_limits<RealT>::denorm_min());

        // The response mode must be an integer inside {0, 1, 2}.
        auto real_mode                      = makeData();
        real_mode.parameters[Params::mode]  = 2.0;
        success                            *= (verifyData(real_mode) > 0);

        auto invalid_mode                      = makeData();
        invalid_mode.parameters[Params::mode]  = static_cast<IdxT>(3);
        success                               *= (verifyData(invalid_mode) > 0);

        // Exercise the serialized integers directly. Named-enum casts here
        // would follow an accidental enum reorder and miss a wire-format bug.
        for (const IdxT mode : {static_cast<IdxT>(0),
                                static_cast<IdxT>(1),
                                static_cast<IdxT>(2)})
        {
          auto mode_data                      = makeData();
          mode_data.parameters[Params::mode]  = mode;
          success                            *= (verifyData(mode_data) == 0);
        }

        // Equal configured limits are valid for every response mode; reversed
        // limits are not.
        for (const IdxT mode : {static_cast<IdxT>(0),
                                static_cast<IdxT>(1),
                                static_cast<IdxT>(2)})
        {
          auto equal                      = makeData();
          equal.parameters[Params::mode]  = mode;
          equal.parameters[Params::Vmin]  = 0.5;
          equal.parameters[Params::Vmax]  = 0.5;
          success                        *= (verifyData(equal) == 0);
        }

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
                                {{index(Vars::XFLOW), 0.8}},
                                "integer-valued component base");

        auto bad_numeric_type                    = makeData();
        bad_numeric_type.parameters[Params::T1]  = true;
        success                                 *= (verifyData(bad_numeric_type) > 0);

        for (const RealT system_base : {-100.0e6,
                                        ZERO<RealT>,
                                        std::numeric_limits<RealT>::denorm_min(),
                                        std::numeric_limits<RealT>::quiet_NaN(),
                                        std::numeric_limits<RealT>::infinity()})
        {
          Fixture<ScalarT> invalid_base(makeData(), system_base);
          success *= (invalid_base.gastpti.allocate() == 0);
          success *= (invalid_base.gastpti.verify() > 0);
        }

        success *= unlinkedSignalRejected<Ext::OMEGA>();
        success *= unlinkedSignalRejected<Ext::PREF>();
        success *= aliasedSignalsRejected();

        Fixture<ScalarT> unallocated(makeData());
        success *= (unallocated.gastpti.initialize() != 0);

        struct TimeConstantCase
        {
          RealT value;
          RealT expected_residual;
        };

        for (const auto& test_case : std::array<TimeConstantCase, 4>{
                 {{0.0, 1.0}, {0.0005, 1.0}, {0.001, 1.0}, {0.002, 0.5}}})
        {
          auto time_data                   = makeData();
          time_data.parameters[Params::T1] = test_case.value;
          time_data.parameters[Params::T2] = test_case.value;
          time_data.parameters[Params::T3] = test_case.value;

          Fixture<ScalarT> time_fixture(time_data);
          success *= time_fixture.initialize(0.4);
          setState(time_fixture.gastpti,
                   {{index(Vars::XVALVE), 0.401},
                    {index(Vars::XFLOW), 0.4},
                    {index(Vars::XTEMP), 0.399},
                    {index(Vars::VLV), 0.402}});
          success *= (time_fixture.evaluate() == 0);
          success *= residualsMatch(
              time_fixture.gastpti,
              {{index(Vars::XVALVE), test_case.expected_residual},
               {index(Vars::XFLOW), test_case.expected_residual},
               {index(Vars::XTEMP), test_case.expected_residual}},
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
        fixture.input(index(Ext::PREF))  = 99.0; // stale value the publication must replace
        success                         *= fixture.initialize(0.4);
        success                         *= (fixture.gastpti.tagDifferentiable() == 0);
        success                         *= (fixture.evaluate() == 0);

        const auto* y  = fixture.gastpti.y().getData();
        success       *= scalarMatches(y[index(Vars::XVALVE)], 0.8, "XVALVE on component base");
        success       *= scalarMatches(y[index(Vars::XFLOW)], 0.8, "XFLOW on component base");
        success       *= scalarMatches(y[index(Vars::XTEMP)], 0.8, "XTEMP on component base");
        success       *= scalarMatches(y[index(Vars::VLOAD)], 0.8, "VLOAD behind the LV gate");
        success       *= scalarMatches(y[index(Vars::VTEMP)], 2.36, "VTEMP at the temperature limit");
        success       *= scalarMatches(y[index(Vars::VLV)], 0.8, "VLV at the fuel flow");
        success       *= scalarExactlyMatches(fixture.pmech(), 0.4, "preserved pmech seed");

        success *= scalarExactlyMatches(fixture.input(index(Ext::OMEGA)),
                                        0.0,
                                        "preserved omega input");
        success *= scalarMatches(fixture.input(index(Ext::PREF)), 0.4, "published pref");

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
          const bool expected = i <= index(Vars::XTEMP);
          if (fixture.gastpti.tag()[i] != expected)
          {
            std::cout << "GASTPTI differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= allResidualsZero(fixture.gastpti);

        // A system-base reference step lands on the droop row scaled by the
        // base ratio.
        fixture.input(index(Ext::PREF))  = 0.5; // the published 0.4 plus a 0.1 step
        success                         *= (fixture.evaluate() == 0);
        success                         *= residualsMatch(fixture.gastpti,
                                                          {{index(Vars::VLOAD), 0.01}},
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
          speed_fixture.input(index(Ext::OMEGA))  = test_case.omega;
          success                                *= speed_fixture.initialize(0.4);
          success                                *= stateMatches(speed_fixture.gastpti,
                                                                 {{index(Vars::XVALVE), test_case.xflow},
                                                                  {index(Vars::XFLOW), test_case.xflow},
                                                                  {index(Vars::XTEMP), test_case.xflow},
                                                                  {index(Vars::VLOAD), test_case.xflow},
                                                                  {index(Vars::VTEMP), test_case.vtemp},
                                                                  {index(Vars::VLV), test_case.xflow}},
                                  "signed nonzero-speed initialization");
          success                                *= scalarExactlyMatches(speed_fixture.pmech(),
                                          0.4,
                                          "signed-speed pmech preservation");
          success                                *= scalarExactlyMatches(speed_fixture.input(index(Ext::OMEGA)),
                                          test_case.omega,
                                          "signed-speed input preservation");
          success                                *= scalarMatches(speed_fixture.input(index(Ext::PREF)),
                                   test_case.pref,
                                   "signed-speed pref publication");
          success                                *= (speed_fixture.evaluate() == 0);
          success                                *= allResidualsZero(speed_fixture.gastpti);
        }

        return success.report(__func__);
      }

      /// Initialization-domain boundaries, response policies, and exact
      /// failure atomicity.
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

        for (const RealT pmech : {std::numeric_limits<RealT>::quiet_NaN(),
                                  std::numeric_limits<RealT>::infinity(),
                                  -std::numeric_limits<RealT>::infinity()})
        {
          success *= initializationRejectedAtomically(
              makeResidualData(), pmech, "non-finite pmech seed");
        }

        for (const RealT omega : {std::numeric_limits<RealT>::quiet_NaN(),
                                  std::numeric_limits<RealT>::infinity(),
                                  -std::numeric_limits<RealT>::infinity()})
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

        // Serialized mode 0 expands its response bounds to an over-rated
        // initial flow and remains exactly at rest.
        auto normal_data                     = makeResidualData();
        normal_data.parameters[Params::mode] = static_cast<IdxT>(0);
        Fixture<ScalarT> over_rated(normal_data);
        over_rated.attachAllInputs();
        success *= over_rated.initialize(0.6); // fuel flow 1.2 above Vmax = 1.1
        success *= stateMatches(over_rated.gastpti,
                                {{index(Vars::XVALVE), 1.2}, {index(Vars::XFLOW), 1.2}, {index(Vars::VLV), 1.2}},
                                "over-rated dispatch");
        success *= scalarExactlyMatches(over_rated.pmech(),
                                        0.6,
                                        "preserved over-rated pmech seed");
        success *= (over_rated.evaluate() == 0);
        success *= allResidualsZero(over_rated.gastpti);

        // Serialized mode 1 derives finite, ordered response bounds with the
        // initial flow below, at, and just above the configured lower limit.
        auto down_only_data                     = makeResidualData();
        down_only_data.parameters[Params::mode] = static_cast<IdxT>(1);
        for (const RealT pmech : {0.0, 0.025, 0.026})
        {
          Fixture<ScalarT> down_only(down_only_data);
          down_only.attachAllInputs();
          success *= down_only.initialize(pmech);
          success *= (down_only.evaluate() == 0);
          success *= allResidualsZero(down_only.gastpti);
        }

        // A failed reinitialization must preserve the last committed response
        // policy as well as state, derivatives, and pref.
        for (const RealT initial_pmech : {0.4, 0.0})
        {
          Fixture<ScalarT> reinitialize(down_only_data);
          reinitialize.attachAllInputs();
          success *= reinitialize.initialize(initial_pmech);
          // Both response policies reach an invalid temperature margin.
          reinitialize.seedPmech(1.0);
          const auto y_before    = copyVector(reinitialize.gastpti.y());
          const auto yp_before   = copyVector(reinitialize.gastpti.yp());
          const auto pref_before = reinitialize.input(index(Ext::PREF));
          if (reinitialize.gastpti.initialize() == 0)
          {
            std::cout << "Expected failed GASTPTI reinitialization\n";
            success = false;
          }
          success *= vectorUnchanged(reinitialize.gastpti.y(), y_before, "reinitialized state");
          success *= vectorUnchanged(reinitialize.gastpti.yp(), yp_before, "reinitialized derivative");
          success *= scalarExactlyMatches(reinitialize.input(index(Ext::PREF)),
                                          pref_before,
                                          "reinitialized pref");

          reinitialize.seedPmech(initial_pmech);
          const RealT xflow = 2.0 * initial_pmech;
          setState(reinitialize.gastpti,
                   {{index(Vars::XVALVE), xflow},
                    {index(Vars::VLV), xflow + 0.25}});
          setDerivative(reinitialize.gastpti, {{index(Vars::XVALVE), 0.0}});
          success        *= (reinitialize.evaluate() == 0);
          RealT expected  = 0.35714285714285715;
          if (initial_pmech == ZERO<RealT>)
          {
            expected = ZERO<RealT>;
          }
          success *= residualsMatch(reinitialize.gastpti,
                                    {{index(Vars::XVALVE), expected}},
                                    "failed reinitialization preserves response policy");
        }

        struct PinnedValveTemperatureCase
        {
          IdxT        mode;
          const char* label;
          RealT       at;
          RealT       vmin;
          RealT       vmax;
          RealT       vtemp;
          RealT       vlv;
        };

        const std::array<PinnedValveTemperatureCase, 4> pinned_temperature_cases{{
            {static_cast<IdxT>(2),
             "Fixed at the temperature limit",
             0.8,
             0.3,
             0.3,
             0.8,
             0.797111886747667},
            {static_cast<IdxT>(2),
             "Fixed above the temperature limit",
             0.0,
             0.3,
             0.3,
             -0.32,
             -0.32000000000000006},
            {static_cast<IdxT>(1),
             "collapsed Down Only at the temperature limit",
             0.8,
             0.8,
             1.1,
             0.8,
             0.797111886747667},
            {static_cast<IdxT>(1),
             "collapsed Down Only above the temperature limit",
             0.0,
             0.8,
             1.1,
             -0.32,
             -0.32000000000000006},
        }};
        for (const auto& test_case : pinned_temperature_cases)
        {
          auto pinned_data                     = makeResidualData();
          pinned_data.parameters[Params::mode] = test_case.mode;
          pinned_data.parameters[Params::At]   = test_case.at;
          pinned_data.parameters[Params::Vmin] = test_case.vmin;
          pinned_data.parameters[Params::Vmax] = test_case.vmax;

          Fixture<ScalarT> pinned(pinned_data);
          pinned.attachAllInputs();
          success *= pinned.initialize(0.4);
          success *= stateMatches(pinned.gastpti,
                                  {{index(Vars::XVALVE), 0.8},
                                   {index(Vars::XFLOW), 0.8},
                                   {index(Vars::XTEMP), 0.8},
                                   {index(Vars::VLOAD), 0.8},
                                   {index(Vars::VTEMP), test_case.vtemp},
                                   {index(Vars::VLV), test_case.vlv},
                                   {index(Vars::PMECH), 0.4}},
                                  test_case.label);
          success *= (pinned.evaluate() == 0);
          success *= allResidualsZero(pinned.gastpti);
        }

        // An unattached reference retains its last successful latch when a
        // later active reinitialization is rejected.
        Fixture<ScalarT> fallback_reinitialize(down_only_data);
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
                                {{index(Vars::XFLOW), 0.0}, {index(Vars::VTEMP), 2.52}},
                                "zero seed");
        success *= (zero_seed.evaluate() == 0);
        success *= allResidualsZero(zero_seed.gastpti);

        auto negative_data                     = makeResidualData();
        negative_data.parameters[Params::Vmin] = -1.0;
        Fixture<ScalarT> negative_seed(negative_data);
        negative_seed.attachAllInputs();
        success *= negative_seed.initialize(-0.1);
        success *= stateMatches(negative_seed.gastpti,
                                {{index(Vars::XFLOW), -0.2},
                                 {index(Vars::VLOAD), -0.2},
                                 {index(Vars::PMECH), -0.1}},
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
                                {{index(Vars::VLOAD), 0.8155903227031184},
                                 {index(Vars::VTEMP), 0.8001000000000001},
                                 {index(Vars::VLV), 0.8}},
                                "near-gate initialization");
        success *= scalarMatches(fixture.input(index(Ext::PREF)),
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
                                {{index(Vars::VLOAD), 0.8},
                                 {index(Vars::VLV), 0.8}},
                                "large finite temperature margin");
        success *= (large_margin.evaluate() == 0);
        success *= allResidualsZero(large_margin.gastpti);

        return success.report(__func__);
      }

      /// Fixed numerical literals provide an independent answer key for all
      /// seven residual rows.
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
        const std::array<Row, index(Vars::MAXIMUM)> expected{{
            {index(Vars::XVALVE), 0.24614285714285705},
            {index(Vars::XFLOW), 0.17755555555555544},
            {index(Vars::XTEMP), -0.001181818181818159},
            {index(Vars::VLOAD), -0.0326},
            {index(Vars::VTEMP), 0.978},
            {index(Vars::VLV), 0.12},
            {index(Vars::PMECH), -0.11239999999999999},
        }};

        success *= (static_cast<size_t>(fixture.gastpti.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.gastpti, expected);

        return success.report(__func__);
      }

      /// Valve anti-windup, speed/damping signs, adjusted Normal limits, and
      /// distinct behavior for serialized modes 0 Normal, 1 Down Only, and
      /// 2 Fixed.
      TestOutcome governorControl()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing GASTPTI response modes and adjusted limits. "
                         "Logged response-limit warnings are expected.");

        // Both response limits block outward motion and admit restoring motion.
        struct AntiWindupCase
        {
          const char* label;
          RealT       xvalve;
          RealT       vlv;
          RealT       expected;
        };

        for (const auto& test_case : std::array<AntiWindupCase, 4>{{
                 {"Vmax blocks an outward valve rate", 1.6, 1.85, 0.0},
                 {"Vmin blocks an outward valve rate", -0.45, -0.7, 0.0},
                 {"Vmax admits a restoring valve rate", 1.6, 1.35, -0.7142857142857143},
                 {"Vmin admits a restoring valve rate", -0.45, -0.2, 0.7142857142857143},
             }})
        {
          Fixture<ScalarT> antiwindup(makeResidualData());
          antiwindup.attachAllInputs();
          success *= antiwindup.initialize(0.4);
          setState(antiwindup.gastpti,
                   {{index(Vars::XVALVE), test_case.xvalve}, {index(Vars::VLV), test_case.vlv}});
          setDerivative(antiwindup.gastpti, {{index(Vars::XVALVE), 0.0}});
          success *= (antiwindup.evaluate() == 0);
          success *= residualsMatch(antiwindup.gastpti,
                                    {{index(Vars::XVALVE), test_case.expected}},
                                    test_case.label);
        }

        // A speed deviation enters the droop and turbine-damping rows.
        Fixture<ScalarT> speed_step(makeResidualData());
        speed_step.attachAllInputs();
        success                             *= speed_step.initialize(0.4);
        speed_step.input(index(Ext::OMEGA))  = 0.05;
        success                             *= (speed_step.evaluate() == 0);
        success                             *= residualsMatch(speed_step.gastpti,
                                                              {{index(Vars::VLOAD), -0.05},
                                                               {index(Vars::PMECH), -0.006}},
                                  "speed deviation in the droop and damping rows");

        struct ModeCase
        {
          IdxT        serialized_mode;
          const char* label;
          RealT       outward_residual;
          RealT       restoring_residual;
        };

        // At the initialized flow xV = 0.8, an upward command distinguishes
        // all three wire values. Down Only is at its derived upper response
        // limit and therefore has the smooth one-half boundary response.
        const std::array<ModeCase, 3> modes{{
            {static_cast<IdxT>(0), "mode 0 Normal", 0.7142857142857143, -0.7142857142857143},
            {static_cast<IdxT>(1), "mode 1 Down Only", 0.35714285714285715, -0.7142857142857143},
            {static_cast<IdxT>(2), "mode 2 Fixed", 0.0, 0.0},
        }};
        for (const auto& test_case : modes)
        {
          auto mode_data                     = makeResidualData();
          mode_data.parameters[Params::mode] = test_case.serialized_mode;
          Fixture<ScalarT> mode_fixture(mode_data);
          mode_fixture.attachAllInputs();
          success *= mode_fixture.initialize(0.4);

          setState(mode_fixture.gastpti,
                   {{index(Vars::XVALVE), 0.8}, {index(Vars::VLV), 1.05}});
          setDerivative(mode_fixture.gastpti, {{index(Vars::XVALVE), 0.0}});
          success *= (mode_fixture.evaluate() == 0);
          success *= residualsMatch(mode_fixture.gastpti,
                                    {{index(Vars::XVALVE), test_case.outward_residual}},
                                    test_case.label);

          setState(mode_fixture.gastpti, {{index(Vars::VLV), 0.55}});
          success *= (mode_fixture.evaluate() == 0);
          success *= residualsMatch(mode_fixture.gastpti,
                                    {{index(Vars::XVALVE), test_case.restoring_residual}},
                                    test_case.label);
        }

        // Pin the smooth Down Only transition just inside, at, and just
        // outside the initialized upper response limit.
        auto down_only_data                     = makeResidualData();
        down_only_data.parameters[Params::mode] = static_cast<IdxT>(1);
        Fixture<ScalarT> down_only(down_only_data);
        down_only.attachAllInputs();
        success                *= down_only.initialize(0.4);
        const RealT transition  = ONE<RealT> / Math::MU<RealT>;
        for (const auto& [offset, expected] : std::array<std::pair<RealT, RealT>, 3>{{
                 {-transition, 0.5221846990214315},
                 {ZERO<RealT>, 0.35714285714285715},
                 {transition, 0.19210101526428275},
             }})
        {
          const RealT xvalve = 0.8 + offset;
          setState(down_only.gastpti,
                   {{index(Vars::XVALVE), xvalve},
                    {index(Vars::VLV), xvalve + 0.25}});
          setDerivative(down_only.gastpti, {{index(Vars::XVALVE), 0.0}});
          success *= (down_only.evaluate() == 0);
          success *= residualsMatch(down_only.gastpti,
                                    {{index(Vars::XVALVE), expected}},
                                    "Down Only upper-response transition");
        }

        // Normal mode expands both sides of the configured interval to admit
        // the initialized flow. The derived boundary must be used thereafter.
        for (const auto& [pmech, command, expected] :
             std::array<std::tuple<RealT, RealT, RealT>, 2>{{
                 {0.6, 0.25, 0.35714285714285715},
                 {0.0, -0.25, -0.35714285714285715},
             }})
        {
          auto normal_data                     = makeResidualData();
          normal_data.parameters[Params::mode] = static_cast<IdxT>(0);
          Fixture<ScalarT> normal(normal_data);
          normal.attachAllInputs();
          success              *= normal.initialize(pmech);
          const RealT boundary  = 2.0 * pmech;
          setState(normal.gastpti,
                   {{index(Vars::XVALVE), boundary},
                    {index(Vars::VLV), boundary + command}});
          setDerivative(normal.gastpti, {{index(Vars::XVALVE), 0.0}});
          success *= (normal.evaluate() == 0);
          success *= residualsMatch(normal.gastpti,
                                    {{index(Vars::XVALVE), expected}},
                                    "Normal adjusted response boundary");
        }

        // Fixed mode pins the valve while the downstream fuel-flow and
        // exhaust-temperature lags remain live.
        auto fixed_data                     = makeResidualData();
        fixed_data.parameters[Params::mode] = static_cast<IdxT>(2);
        Fixture<ScalarT> fixed(fixed_data);
        fixed.attachAllInputs();
        success *= fixed.initialize(0.4);
        setAnswerKeyInputs(fixed);
        setAnswerKeyState(fixed.gastpti);
        success *= (fixed.evaluate() == 0);
        success *= residualsMatch(fixed.gastpti,
                                  {{index(Vars::XVALVE), -0.011},
                                   {index(Vars::XFLOW), 0.17755555555555544},
                                   {index(Vars::XTEMP), -0.001181818181818159}},
                                  "mode 2 pins only the fuel valve");

        // From initialized rest, a speed step still enters both live algebraic
        // interface rows.
        Fixture<ScalarT> fixed_interface(fixed_data);
        fixed_interface.attachAllInputs();
        success                                  *= fixed_interface.initialize(0.4);
        fixed_interface.input(index(Ext::OMEGA))  = 0.05;
        success                                  *= (fixed_interface.evaluate() == 0);
        success                                  *= residualsMatch(
            fixed_interface.gastpti,
            {{index(Vars::XVALVE), 0.0},
                                              {index(Vars::XFLOW), 0.0},
                                              {index(Vars::XTEMP), 0.0},
                                              {index(Vars::VLOAD), -0.05},
                                              {index(Vars::PMECH), -0.006}},
            "Fixed mode preserves its live algebraic interface");

        // Down Only pins the valve when initialization is below or exactly at
        // Vmin. A point just above Vmin retains a narrow active interval where
        // the smooth lower- and upper-limit gates interact at the initialized
        // upper boundary.
        for (const auto& [pmech, expected] :
             std::array<std::pair<RealT, RealT>, 3>{{
                 {0.0, 0.0},
                 {0.025, 0.0},
                 // The 0.002-wide interval is narrower than the 1/MU smooth
                 // transition. At its upper boundary the lower gate is still
                 // 0.61774787476924897, so the combined indicator is
                 // 1 - 0.5 * lower_gate = 0.69112606261537546.
                 {0.026, 0.4936614732966968},
             }})
        {
          Fixture<ScalarT> narrow_down_only(down_only_data);
          narrow_down_only.attachAllInputs();
          success              *= narrow_down_only.initialize(pmech);
          const RealT boundary  = 2.0 * pmech;
          setState(narrow_down_only.gastpti,
                   {{index(Vars::XVALVE), boundary},
                    {index(Vars::VLV), boundary + 0.25}});
          setDerivative(narrow_down_only.gastpti,
                        {{index(Vars::XVALVE), 0.0}});
          success *= (narrow_down_only.evaluate() == 0);
          success *= residualsMatch(
              narrow_down_only.gastpti,
              {{index(Vars::XVALVE), expected}},
              "Down Only response interval at the configured lower limit");
        }

        // A collapsed Down Only interval pins only the valve; downstream lags
        // retain the same restoring equations as Normal mode.
        Fixture<ScalarT> collapsed(down_only_data);
        collapsed.attachAllInputs();
        success *= collapsed.initialize(0.0);
        setAnswerKeyState(collapsed.gastpti);
        success *= (collapsed.evaluate() == 0);
        success *= residualsMatch(collapsed.gastpti,
                                  {{index(Vars::XVALVE), -0.011},
                                   {index(Vars::XFLOW), 0.17755555555555544},
                                   {index(Vars::XTEMP), -0.001181818181818159}},
                                  "collapsed Down Only pins only the fuel valve");

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

        for (const auto& test_case : std::array<GateCase, 3>{{
                 {"the load demand wins the LV gate", 0.3, 1.5, 0.3},
                 {"the temperature demand wins the LV gate", 1.5, 0.3, 0.30000000000000004},
                 {"equal demands split the smooth LV gate", 0.9, 0.9, 0.897111886747667},
             }})
        {
          Fixture<ScalarT> gate(makeResidualData());
          gate.attachAllInputs();
          success *= gate.initialize(0.4);
          setState(gate.gastpti,
                   {{index(Vars::VLOAD), test_case.vload},
                    {index(Vars::VTEMP), test_case.vtemp},
                    {index(Vars::VLV), 0.0}});
          success *= (gate.evaluate() == 0);
          success *= residualsMatch(gate.gastpti,
                                    {{index(Vars::VLV), test_case.expected}},
                                    test_case.label);
        }

        // The exhaust-temperature feedback drives the temperature demand.
        Fixture<ScalarT> feedback(makeResidualData());
        feedback.attachAllInputs();
        success *= feedback.initialize(0.4);
        setState(feedback.gastpti,
                 {{index(Vars::XTEMP), 0.9}, {index(Vars::VTEMP), 1.1}});
        success *= (feedback.evaluate() == 0);
        success *= residualsMatch(feedback.gastpti,
                                  {{index(Vars::VTEMP), 1.06}},
                                  "temperature feedback");

        return success.report(__func__);
      }

      /// Fixed derivative oracles cover all rows and response regions, with
      /// Enzyme agreement checked when available.
      TestOutcome jacobian()
      {
        using DepMap = DependencyTracking::Variable::DependencyMap;
        using DepVar = DependencyTracking::Variable;

        TestStatus success = true;

        const auto data = makeResidualData();

        const auto                interior = dependencyTrackingJacobian(data, success);
        const std::vector<DepMap> interior_expected{
            {{index(Vars::XVALVE), -3.8571428571428572},
             {index(Vars::VLV), 2.8571428571428572}},
            {{index(Vars::XVALVE), 2.2222222222222223},
             {index(Vars::XFLOW), -3.2222222222222223}},
            {{index(Vars::XFLOW), 0.45454545454545453},
             {index(Vars::XTEMP), -1.4545454545454546}},
            {{index(Vars::VLOAD), -0.06},
             {index(Vars::MAXIMUM) + index(Ext::OMEGA), -1.0},
             {index(Vars::MAXIMUM) + index(Ext::PREF), 0.12}},
            {{index(Vars::XTEMP), -0.4},
             {index(Vars::VTEMP), -1.0}},
            {{index(Vars::VLOAD), 1.0},
             {index(Vars::VLV), -1.0}},
            {{index(Vars::XFLOW), 1.0},
             {index(Vars::PMECH), -2.0},
             {index(Vars::MAXIMUM) + index(Ext::OMEGA), -0.12}},
        };
        success *= jacobianMatches(interior,
                                   interior_expected,
                                   "Normal interior answer key");

        // At equality, the smooth low-value selector splits its sensitivity
        // evenly between the two demand signals.
        {
          Fixture<DepVar> selector(data);
          selector.attachAllInputs();
          success *= selector.initialize(0.4);
          setState(selector.gastpti,
                   {{index(Vars::VLOAD), 0.9},
                    {index(Vars::VTEMP), 0.9},
                    {index(Vars::VLV), 0.7}});
          numberVariables(selector);
          success *= (selector.evaluate() == 0);

          const DepMap expected{{index(Vars::VLOAD), 0.5},
                                {index(Vars::VTEMP), 0.5},
                                {index(Vars::VLV), -1.0}};
          success *= jacobianRowMatches(
              selector.gastpti.getResidual().getData()[index(Vars::VLV)].getDependencies(),
              expected,
              index(Vars::VLV),
              "selector equality");
        }

        // The one-half response at an active upper bound has a large state
        // derivative from the smooth gate. Exercise both a Normal interval
        // expanded by initialization and the Down Only upper boundary.
        success *= responseBoundaryJacobian(data,
                                            static_cast<IdxT>(0),
                                            0.6,
                                            1.2,
                                            "Normal adjusted upper bound");
        success *= responseBoundaryJacobian(data,
                                            static_cast<IdxT>(1),
                                            0.4,
                                            0.8,
                                            "Down Only upper bound");

        // Fixed mode removes only the valve-drive coupling. The downstream
        // lag rows and all four algebraic rows remain unchanged.
        auto fixed_data                     = data;
        fixed_data.parameters[Params::mode] = static_cast<IdxT>(2);
        const auto                fixed     = dependencyTrackingJacobian(fixed_data, success);
        const std::vector<DepMap> fixed_expected{
            {{index(Vars::XVALVE), -1.0}},
            interior_expected[index(Vars::XFLOW)],
            interior_expected[index(Vars::XTEMP)],
            interior_expected[index(Vars::VLOAD)],
            interior_expected[index(Vars::VTEMP)],
            interior_expected[index(Vars::VLV)],
            interior_expected[index(Vars::PMECH)],
        };
        success *= jacobianMatches(fixed,
                                   fixed_expected,
                                   "Fixed-mode answer key");

#ifdef GRIDKIT_ENABLE_ENZYME
        const auto enzyme  = enzymeJacobian(data, success);
        success           *= jacobianMatches(enzyme,
                                   interior,
                                   "Enzyme versus dependency tracking");
#endif

        return success.report(__func__);
      }

    private:
      using Params = PhasorDynamics::Governor::GastPtiParameters;
      using Vars   = PhasorDynamics::Governor::GastPtiInternalVariables;
      using Ext    = PhasorDynamics::Governor::GastPtiExternalVariables;
      using Mon    = PhasorDynamics::Governor::GastPtiMonitorableVariables;
      using Mode   = PhasorDynamics::Governor::ResponseMode;
      using Data   = PhasorDynamics::Governor::GastPtiData<RealT, IdxT>;

      static constexpr size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(Ext variable)
      {
        return static_cast<size_t>(variable);
      }

      /// A vector row paired with a value: either an input to write or an
      /// expected result. Rows are derived from the canonical variable enums.
      using Row      = std::pair<size_t, RealT>;
      using Rows     = std::initializer_list<Row>;
      using GastPtiT = PhasorDynamics::Governor::GastPti<ScalarT, IdxT>;

      /// Owns the GASTPTI model, the assigned mechanical-power node, and the
      /// attached input nodes. Signal storage is declared before the model so
      /// every referenced node outlives GASTPTI. Copying would invalidate the
      /// model and signal-node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, index(Ext::MAXIMUM)>                                   input_values_{};
        std::array<IdxT, index(Ext::MAXIMUM)>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, index(Ext::MAXIMUM)> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> pmech_node_;

      public:
        explicit Fixture(const Data& data, RealT system_va_base = 100.0e6)
          : gastpti(data)
        {
          gastpti.setSystemBase(60.0, system_va_base);
          gastpti.getSignals().template assignSignalNode<Vars::PMECH>(&pmech_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Attach fixture-owned storage to every external input.
        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = gastpti.size();

          for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = gastpti.getSignals();
          signals.template attachSignalNode<Ext::OMEGA>(&input_nodes_[index(Ext::OMEGA)]);
          signals.template attachSignalNode<Ext::PREF>(&input_nodes_[index(Ext::PREF)]);
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
        data.parameters[Params::mode]  = static_cast<IdxT>(Mode::Normal);
        return data;
      }

      Data makeData() const
      {
        auto data = makeMinimalData();

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
        data.parameters[Params::mode]  = static_cast<IdxT>(Mode::Normal);
        return data;
      }

      Data makeResidualData() const
      {
        auto data = makeData();

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
        fixture.input(index(Ext::OMEGA)) = static_cast<T>(0.02);
        fixture.input(index(Ext::PREF))  = static_cast<T>(0.31);
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti) const
      {
        setState(gastpti,
                 {{index(Vars::XVALVE), 0.62},
                  {index(Vars::XFLOW), 0.55},
                  {index(Vars::XTEMP), 0.48},
                  {index(Vars::VLOAD), 0.83},
                  {index(Vars::VTEMP), 1.35},
                  {index(Vars::VLV), 0.71},
                  {index(Vars::PMECH), 0.33}});
        setDerivative(gastpti,
                      {{index(Vars::XVALVE), 0.011},
                       {index(Vars::XFLOW), -0.022},
                       {index(Vars::XTEMP), 0.033}});
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

      bool invalidParameterCase(Params parameter, RealT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        return verifyData(data) > 0;
      }

      int verifyData(const Data& data) const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>        pmech;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> model(data);
        model.getSignals().template assignSignalNode<Vars::PMECH>(&pmech);
        return model.verify();
      }

      template <Ext variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>        unlinked_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>        pmech_node;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> model(makeData());
        model.getSignals().template assignSignalNode<Vars::PMECH>(&pmech_node);
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
        pref_alias.getSignals().template assignSignalNode<Vars::PMECH>(&pmech_pref);
        pref_alias.getSignals().template attachSignalNode<Ext::PREF>(&pmech_pref);
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
        speed_alias.getSignals().template assignSignalNode<Vars::PMECH>(&pmech_speed);
        speed_alias.getSignals().template attachSignalNode<Ext::OMEGA>(&pmech_speed);
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
        input_alias.getSignals().template assignSignalNode<Vars::PMECH>(&pmech);
        input_alias.getSignals().template attachSignalNode<Ext::OMEGA>(&shared_input);
        input_alias.getSignals().template attachSignalNode<Ext::PREF>(&shared_input);
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
          if (!exactlyEqual(static_cast<RealT>(values[i]), snapshot[i]))
          {
            std::cout << "GASTPTI " << what << " row " << i
                      << " changed: "
                      << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                      << static_cast<RealT>(values[i]) << " != " << snapshot[i]
                      << '\n';
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
        fixture.input(index(Ext::OMEGA)) = omega;
        fixture.input(index(Ext::PREF))  = 77.0; // must stay untouched on rejection
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

        if (!scalarExactlyMatches(fixture.pmech(),
                                  pmech,
                                  "rejected pmech seed preservation"))
        {
          success = false;
        }
        if (!scalarExactlyMatches(fixture.input(index(Ext::OMEGA)),
                                  omega,
                                  "rejected omega preservation"))
        {
          success = false;
        }
        if (!scalarExactlyMatches(fixture.input(index(Ext::PREF)),
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
      void setState(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti, Rows rows) const
      {
        auto* y = gastpti.y().getData();
        for (const auto& [row, value] : rows)
        {
          y[row] = static_cast<T>(value);
        }
        gastpti.y().setDataUpdated();
      }

      /// setState() for the derivative vector.
      template <typename T>
      void setDerivative(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti, Rows rows) const
      {
        auto* yp = gastpti.yp().getData();
        for (const auto& [row, value] : rows)
        {
          yp[row] = static_cast<T>(value);
        }
        gastpti.yp().setDataUpdated();
      }

      /// Compare one vector row against its expected value. Every row check
      /// in this suite reports through here, so failures share one format.
      /// Rows are named by the canonical enum position used by the expectation.
      static bool rowMatches(RealT       actual,
                             RealT       expected,
                             const char* what,
                             size_t      row,
                             const char* context)
      {
        if (isEqual(actual, expected, kBehaviorTol))
        {
          return true;
        }
        std::cout << "GASTPTI " << what << " row " << row << ' ' << context
                  << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual
                  << " != " << expected << '\n';
        return false;
      }

      /// Check selected rows of a model vector against expected values.
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
          if (!rowMatches(static_cast<RealT>(values[row]), expected, what, row, context))
          {
            success = false;
          }
        }
        return success;
      }

      bool residualsMatch(const GastPtiT& gastpti, Rows rows, const char* context = "") const
      {
        return rowsMatch(gastpti.getResidual(), rows.begin(), rows.size(), "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const GastPtiT&              gastpti,
                          const std::array<Row, size>& rows,
                          const char*                  context = "") const
      {
        return rowsMatch(gastpti.getResidual(), rows.data(), size, "residual", context);
      }

      bool stateMatches(const GastPtiT& gastpti, Rows rows, const char* context = "") const
      {
        return rowsMatch(gastpti.y(), rows.begin(), rows.size(), "state", context);
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
                         ScalarT     tolerance = kBehaviorTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual
                  << " != " << expected << "\n";
        return false;
      }

      static bool exactlyEqual(RealT actual, RealT expected)
      {
        return std::bit_cast<std::array<std::byte, sizeof(RealT)>>(actual)
               == std::bit_cast<std::array<std::byte, sizeof(RealT)>>(expected);
      }

      bool scalarExactlyMatches(ScalarT     actual,
                                ScalarT     expected,
                                const char* label) const
      {
        const RealT actual_value   = static_cast<RealT>(actual);
        const RealT expected_value = static_cast<RealT>(expected);
        if (exactlyEqual(actual_value, expected_value))
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

      using DependencyMap = DependencyTracking::Variable::DependencyMap;

      static RealT dependencyValue(const DependencyMap& row, size_t column)
      {
        const auto entry = row.find(column);
        if (entry == row.end())
        {
          return ZERO<RealT>;
        }
        return entry->second;
      }

      bool responseBoundaryJacobian(const Data& data,
                                    IdxT        serialized_mode,
                                    RealT       pmech,
                                    RealT       boundary,
                                    const char* label) const
      {
        using DepVar = DependencyTracking::Variable;

        auto boundary_data                     = data;
        boundary_data.parameters[Params::mode] = serialized_mode;
        Fixture<DepVar> fixture(boundary_data);
        fixture.attachAllInputs();

        bool success = fixture.initialize(pmech);
        setState(fixture.gastpti,
                 {{index(Vars::XVALVE), boundary},
                  {index(Vars::VLV), boundary + 0.25}});
        setDerivative(fixture.gastpti, {{index(Vars::XVALVE), 0.0}});
        numberVariables(fixture);
        if (fixture.evaluate() != 0)
        {
          success = false;
        }

        const DependencyMap expected{
            {index(Vars::XVALVE), -45.285714285714285},
            {index(Vars::VLV), 1.4285714285714286},
        };
        if (!jacobianRowMatches(
                fixture.gastpti.getResidual().getData()[index(Vars::XVALVE)].getDependencies(),
                expected,
                index(Vars::XVALVE),
                label))
        {
          success = false;
        }
        return success;
      }

      bool jacobianRowMatches(const DependencyMap& actual,
                              const DependencyMap& expected,
                              size_t               row,
                              const char*          context) const
      {
        bool success = true;

        for (const auto& [column, value] : actual)
        {
          static_cast<void>(value);
          if (!jacobianColumnMatches(actual, expected, row, column, context))
          {
            success = false;
          }
        }
        for (const auto& [column, value] : expected)
        {
          static_cast<void>(value);
          if (!actual.contains(column)
              && !jacobianColumnMatches(actual, expected, row, column, context))
          {
            success = false;
          }
        }
        return success;
      }

      bool jacobianColumnMatches(const DependencyMap& actual,
                                 const DependencyMap& expected,
                                 size_t               row,
                                 size_t               column,
                                 const char*          context) const
      {
        const RealT actual_value   = dependencyValue(actual, column);
        const RealT expected_value = dependencyValue(expected, column);
        if (isEqual(actual_value, expected_value, kJacobianTol))
        {
          return true;
        }

        std::cout << "GASTPTI Jacobian row " << row << ", column "
                  << column << " " << context << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual_value << " != " << expected_value << '\n';
        return false;
      }

      bool jacobianMatches(const std::vector<DependencyMap>& actual,
                           const std::vector<DependencyMap>& expected,
                           const char*                       context) const
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
          if (!jacobianRowMatches(actual[row], expected[row], row, context))
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
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
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

        const auto                         model_size = static_cast<size_t>(fixture.gastpti.size());
        std::vector<DepVar::DependencyMap> rows(model_size);
        const auto*                        f = fixture.gastpti.getResidual().getData();
        for (size_t i = 0; i < model_size; ++i)
        {
          rows[i] = f[i].getDependencies();
        }
        return rows;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
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
