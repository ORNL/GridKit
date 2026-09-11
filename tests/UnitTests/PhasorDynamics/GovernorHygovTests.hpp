#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "ComponentTestFixture.hpp"

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    class GovernorHygovTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      GovernorHygovTests()  = default;
      ~GovernorHygovTests() = default;

      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

      /// Construction and every verify() error class, including parameter
      /// types and finiteness, parameter relationships, power bases, curve
      /// shape, gate-limit domain, the required pmech assignment, and signal
      /// linkage, plus differentiability tagging.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> empty;
        success *= (empty.size() == static_cast<IdxT>(Internal::MAXIMUM));
        success *= (empty.getMonitor() == nullptr);

        Fixture<ScalarT> configured(makeData(), __func__, kTol);
        configureHygov(configured);
        success *= (configured.model().size() == static_cast<IdxT>(Internal::MAXIMUM));
        success *= (configured.model().getMonitor() != nullptr);
        success *= (configured.model().verify() == 0);

        const RestoreVerbosity restore_verbosity;
        // Suppress expected errors and warnings from the invalid cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        Fixture<ScalarT> minimal(makeMinimalData(), __func__, kTol);
        configureHygov(minimal);
        success *= (minimal.model().verify() == 0);
        success *= defaultsMatchDocumentedValues();

        auto missing_trate_data = makeMinimalData();
        missing_trate_data.parameters.erase(Params::Trate);
        Fixture<ScalarT> missing_trate(missing_trate_data, __func__, kTol);
        configureHygov(missing_trate);
        success *= (missing_trate.model().verify() > 0);

        success *= (empty.verify() > 0);

        const RealT nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT infinity = std::numeric_limits<RealT>::infinity();

        const std::array<Params, 30> real_parameters{{
            Params::Trate,
            Params::Rperm,
            Params::Rtemp,
            Params::Tr,
            Params::Tf,
            Params::Tg,
            Params::Velm,
            Params::Gmax,
            Params::Gmin,
            Params::Tw,
            Params::At,
            Params::Dturb,
            Params::Qnl,
            Params::Tn,
            Params::Tnp,
            Params::db1,
            Params::db2,
            Params::Hdam,
            Params::Gv0,
            Params::Gv1,
            Params::Gv2,
            Params::Gv3,
            Params::Gv4,
            Params::Gv5,
            Params::Pgv0,
            Params::Pgv1,
            Params::Pgv2,
            Params::Pgv3,
            Params::Pgv4,
            Params::Pgv5,
        }};
        const std::array<RealT, 3>   nonfinite_values{{nan, infinity, -infinity}};

        for (const Params parameter : real_parameters)
        {
          for (const RealT value : nonfinite_values)
          {
            Fixture<ScalarT> invalid(withParameters(makeData(), {{parameter, value}}), __func__, kTol);
            configureHygov(invalid);
            success *= (invalid.model().verify() > 0);
          }
        }

        // The pmech output is required, so a model without an assigned node
        // is rejected even when every parameter is valid.
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> unassigned(makeData());
        success *= (unassigned.verify() > 0);

        const std::array<std::pair<Params, RealT>, 19> invalid_parameter_values{{
            {Params::Trate, 0.0},
            {Params::Trate, -1.0},
            {Params::Rtemp, 0.0},
            {Params::Tr, -0.1},
            {Params::Tf, -0.1},
            {Params::Tg, -0.1},
            {Params::Tw, -0.1},
            {Params::Tn, -0.1},
            {Params::Tnp, -0.1},
            {Params::Velm, -0.1},
            {Params::Gmin, 1.1},
            {Params::At, 0.0},
            {Params::Dturb, -0.1},
            {Params::db1, -0.1},
            {Params::Hdam, 0.0},
            {Params::Gv2, 0.1},
            {Params::Pgv2, 0.1},
            {Params::Gmin, -0.05},
            {Params::Gmax, 1.05},
        }};

        for (const auto& [parameter, value] : invalid_parameter_values)
        {
          Fixture<ScalarT> invalid(withParameters(makeData(), {{parameter, value}}), __func__, kTol);
          configureHygov(invalid);
          success *= (invalid.model().verify() > 0);
        }

        // A curve with no rise cannot yield a unique gate.
        Fixture<ScalarT> flat_curve(withParameters(makeData(), {{Params::Pgv1, 0.0}, {Params::Pgv2, 0.0}, {Params::Pgv3, 0.0}, {Params::Pgv4, 0.0}, {Params::Pgv5, 0.0}}), __func__, kTol);
        configureHygov(flat_curve);
        success *= (flat_curve.model().verify() > 0);

        // A curve that rises only outside the configured response limits is
        // valid because initialization may expand those limits.
        Fixture<ScalarT> flat_range(withParameters(makeData(), {{Params::Gmin, 0.0}, {Params::Gmax, 0.2}, {Params::Pgv0, 0.5}, {Params::Pgv1, 0.5}, {Params::Pgv2, 0.5}, {Params::Pgv3, 0.5}, {Params::Pgv4, 0.5}, {Params::Pgv5, 1.0}}), __func__, kTol);
        configureHygov(flat_range);
        success *= (flat_range.model().verify() == 0);

        // A requested backlash is accepted, warns, and remains inactive.
        Fixture<ScalarT> backlash(withParameters(makeData(), {{Params::db2, 0.5}}), __func__, kTol);
        configureHygov(backlash);
        success *= (backlash.model().verify() == 0);

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                   = makeData();
        integer_real.parameters[Params::Tw] = static_cast<IdxT>(2);
        Fixture<ScalarT> integer(integer_real, __func__, kTol);
        configureHygov(integer);
        success *= (integer.model().verify() == 0);

        auto bad_numeric_type                      = makeData();
        bad_numeric_type.parameters[Params::Trate] = true;
        Fixture<ScalarT> bad_type(bad_numeric_type, __func__, kTol);
        configureHygov(bad_type);
        success *= (bad_type.model().verify() > 0);

        Fixture<ScalarT> base_overflow(withParameters(makeData(), {{Params::Trate, std::numeric_limits<RealT>::max()}}), __func__, kTol);
        configureHygov(base_overflow);
        success *= (base_overflow.model().verify() > 0);

        Fixture<ScalarT> ratio_overflow(withParameters(makeData(), {{Params::Trate, std::numeric_limits<RealT>::min()}}), __func__, kTol);
        configureHygov(ratio_overflow);
        success *= (ratio_overflow.model().verify() > 0);

        const std::array<RealT, 6> invalid_system_bases{{
            0.0,
            -1.0,
            nan,
            infinity,
            -infinity,
            std::numeric_limits<RealT>::min(),
        }};

        for (const RealT system_base : invalid_system_bases)
        {
          Fixture<ScalarT> invalid_base(makeData(), __func__, kTol);
          configureHygov(invalid_base, system_base);
          success *= (invalid_base.model().verify() > 0);
        }

        success *= unlinkedSignalRejected<External::OMEGA>();
        success *= unlinkedSignalRejected<External::PREF>();
        success *= unlinkedSignalRejected<External::PAUX>();

        // All five zero time constants use the documented numerical floor and
        // still admit a consistent steady-state initialization.
        Fixture<ScalarT> floors(withParameters(makeData(), {{Params::Tr, 0.0}, {Params::Tf, 0.0}, {Params::Tg, 0.0}, {Params::Tw, 0.0}, {Params::Tnp, 0.0}}), __func__, kTol);
        configureHygov(floors);
        if (!floors.initialize({{Internal::PMECH, 0.4}}))
          return TestStatus(false).report(__func__);
        success *= floors.checkSteadyState();

        return success.report(__func__);
      }

      /// A nonidentity power-base initialization with every port attached.
      /// The machine-provided pmech value must remain unchanged while HYGOV
      /// initializes and publishes its resolved load reference.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(withParameters(makeData(), {{Params::Trate, 50.0}}), __func__, kTol);
        configureHygov(fixture);
        attachInputs(fixture);
        if (!fixture.setInput(External::PAUX, 0.02))
          return TestStatus(false).report(__func__);
        if (!fixture.setInput(External::PREF, 99.0))
          return TestStatus(false).report(__func__); // stale value the publication must replace
        if (!fixture.initialize({{Internal::PMECH, 0.4}}))
          return TestStatus(false).report(__func__);
        success *= (fixture.model().tagDifferentiable() == 0);
        if (!fixture.evaluateResidual())
          return TestStatus(false).report(__func__);

        const auto* y  = fixture.model().y().getData();
        success       *= scalarMatches(y[static_cast<size_t>(Internal::XF)], 0.0, "XF at rest");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::C)],
                                 0.9,
                                 "C on component base",
                                 kTol + 2.0 * rampError(0.1));
        success       *= scalarMatches(y[static_cast<size_t>(Internal::G)],
                                 0.9,
                                 "G on component base",
                                 kTol + 2.0 * rampError(0.1));
        success       *= scalarMatches(y[static_cast<size_t>(Internal::Q)], 0.9, "Q on component base");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::PGV)],
                                 0.9,
                                 "PGV on component base");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::H)], 1.0, "H at the dam head");
        success       *= scalarMatches(fixture.output(Internal::PMECH), 0.4, "preserved pmech value");

        success *= scalarMatches(fixture.input(External::OMEGA), 0.0, "preserved omega input");
        success *= scalarMatches(fixture.input(External::PREF), 0.0025, "published pref", kTol + 0.05 * rampError(0.1));
        success *= scalarMatches(fixture.input(External::PAUX), 0.02, "preserved paux input");

        // Verify the six documented outputs through the public monitor controller.
        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.model().getMonitor());
        std::stringstream monitor_output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, monitor_output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string monitor_header;
        std::string monitor_values;
        std::getline(monitor_output, monitor_header);
        std::getline(monitor_output, monitor_values);
        success              *= (monitor_header == "t,Hygov_hygov_test_pmech,Hygov_hygov_test_filter,"
                                                   "Hygov_hygov_test_desiredgate,Hygov_hygov_test_gate,"
                                                   "Hygov_hygov_test_flow,Hygov_hygov_test_head");
        const auto monitored  = Tokenizer<RealT>(monitor_values, ',')();
        if (monitored.size() == 7)
        {
          success *= scalarMatches(monitored[1], 0.4, "monitored pmech");
          success *= scalarMatches(monitored[2], 0.0, "monitored filter");
          success *= scalarMatches(monitored[3], fixture.state(Internal::C), "monitored desiredgate");
          success *= scalarMatches(monitored[4], fixture.state(Internal::G), "monitored gate");
          success *= scalarMatches(monitored[5], 0.9, "monitored flow");
          success *= scalarMatches(monitored[6], 1.0, "monitored head");
        }
        else
        {
          std::cout << "HYGOV monitor emitted " << monitored.size()
                    << " values instead of 7\n";
          success = false;
        }

        // The five governor states carry derivatives; the rest is algebraic.
        for (size_t i = 0; i < static_cast<size_t>(fixture.model().size()); ++i)
        {
          const bool differential = i <= static_cast<size_t>(Internal::Q);
          if (fixture.model().tag()[i] != differential)
          {
            std::cout << "HYGOV differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }

        success *= fixture.checkSteadyState();

        // A system-base reference step lands on the governor error scaled by
        // the base ratio.
        if (!fixture.setInput(External::PREF, fixture.input(External::PREF) + 0.1))
          return TestStatus(false).report(__func__);
        success *= fixture.checkResidualRows({{Internal::EF, 0.2}}, "reference step on the component base");

        // Unattached ports fall back to the references latched by
        // initialize(), so the same steady state holds without a controller.
        Fixture<ScalarT> latched(withParameters(makeData(), {{Params::Trate, 50.0}}), __func__, kTol);
        configureHygov(latched);
        if (!latched.initialize({{Internal::PMECH, 0.4}}))
          return TestStatus(false).report(__func__);
        success *= latched.checkSteadyState();

        return success.report(__func__);
      }

      /// Mechanical-power, gate-limit, speed-deviation, and finite-input
      /// initialization domains. Response limits and high-power dam head
      /// are adjusted when needed; every rejected initialization is atomic.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        const RestoreVerbosity restore_verbosity;
        // Suppress expected errors and limit-adjustment warnings from the cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        success *= initializationRejectedAtomically(
            makeResidualData(),
            -0.3,
            {{External::OMEGA, 0.0}, {External::PREF, 77.0}, {External::PAUX, 0.02}},
            "mechanical power below the gate curve");

        const auto no_finite_head = withParameters(
            makeData(),
            {{Params::Pgv0, -1.0},
             {Params::Pgv1, -0.8},
             {Params::Pgv2, -0.6},
             {Params::Pgv3, -0.4},
             {Params::Pgv4, -0.2},
             {Params::Pgv5, 0.0}});
        success *= initializationRejectedAtomically(
            no_finite_head,
            0.0,
            {{External::OMEGA, 0.0}, {External::PREF, 77.0}, {External::PAUX, 0.02}},
            "no finite effective Hdam");

        // At = 1 and Qnl = 0 give P = H^(3/2): 8 pu needs Q = 2 and H = 4.
        Fixture<ScalarT> adjusted(withParameters(makeData(), {{Params::Trate, 1.8}, {Params::At, 1.0}, {Params::Qnl, 0.0}, {Params::Gmax, 0.5}}), __func__, kTol);
        if (!initializeHygov(adjusted, 0.144))
          return TestStatus(false).report(__func__);
        success          *= adjusted.checkStateRows({{Internal::C, 1.0}, {Internal::G, 1.0}}, "fully open gate");
        const RealT knee  = std::log(2.0) / Math::MU<RealT>;
        success          *= scalarMatches(adjusted.state(Internal::PGV), 1.0, "gate power", kTol + knee);
        success          *= scalarMatches(adjusted.state(Internal::Q), 2.0, "high-power flow", kTol + 2.0 * knee);
        success          *= scalarMatches(adjusted.state(Internal::H), 4.0, "effective dam head", kTol + 2.0 * knee);
        success          *= scalarMatches(adjusted.output(Internal::PMECH), 0.144, "preserved pmech value");
        success          *= scalarMatches(adjusted.input(External::PREF),
                                 0.0009,
                                 "published pref");
        success          *= adjusted.checkSteadyState();

        struct ResponseLimitCase
        {
          const char* label;
          Params      limit_parameter;
          RealT       limit;
          RealT       rate;
        };

        const std::array<ResponseLimitCase, 2> response_limit_cases{{
            {"expanded upper response limit", Params::Gmax, 0.5, 0.1},
            {"expanded lower response limit", Params::Gmin, 0.7, -0.1},
        }};

        for (const auto& test_case : response_limit_cases)
        {
          Fixture<ScalarT> fixture(withParameters(makeResidualData(), {{test_case.limit_parameter, test_case.limit}}), __func__, kTol);
          configureHygov(fixture);
          if (!fixture.initialize({{Internal::PMECH, 0.4}}))
            return TestStatus(false).report(__func__);
          const RealT gate            = fixture.state(Internal::C);
          const bool  gate_is_outside = test_case.rate > 0.0
                                            ? gate > test_case.limit
                                            : gate < test_case.limit;
          if (!gate_is_outside)
          {
            std::cout << test_case.label << " did not initialize outside the configured limit\n";
            success = false;
          }
          success *= fixture.checkStateRows({{Internal::G, gate}, {Internal::H, 1.2}}, test_case.label);
          success *= fixture.checkSteadyState();

          // The effective response bound admits an outward rate between the
          // configured limit and initialized gate.
          if (!fixture.setPoint({.state      = {{Internal::C, 0.5 * (test_case.limit + gate)},
                                                {Internal::RC, test_case.rate}},
                                 .derivative = {{Internal::C, 0.0}}})
              || !fixture.evaluateResidual())
            return TestStatus(false).report(__func__);
          const RealT response_rate    = fixture.residual(Internal::C);
          const bool  rate_is_admitted = test_case.rate > 0.0
                                             ? response_rate > 0.5 * test_case.rate
                                             : response_rate < 0.5 * test_case.rate;
          if (!rate_is_admitted)
          {
            std::cout << test_case.label << " did not admit the outward desired-gate rate\n";
            success = false;
          }
        }

        // A failed retry preserves the effective head and response bounds from
        // the prior success.
        if (!adjusted.setInput(External::OMEGA, 0.03))
          return TestStatus(false).report(__func__);
        const auto effective_before  = adjusted.snapshot();
        success                     *= adjusted.model().initialize() != 0;
        success                     *= adjusted.checkUnchanged(effective_before);
        if (!adjusted.setInput(External::OMEGA, 0.0))
          return TestStatus(false).report(__func__);
        success *= adjusted.checkSteadyState();

        if (!adjusted.setPoint({.state = {{Internal::C, 0.75}, {Internal::RC, 0.2}}, .derivative = {{Internal::C, 0.0}}})
            || !adjusted.evaluateResidual())
          return TestStatus(false).report(__func__);
        const RealT preserved_rate = adjusted.residual(Internal::C);
        if (!(preserved_rate > 0.19))
        {
          std::cout << "failed initialization did not preserve effective Gmax\n";
          success = false;
        }

        // A later feasible initialization starts again from configured limits
        // and Hdam.
        if (!adjusted.setState({{Internal::PMECH, 0.009}}))
          return TestStatus(false).report(__func__);
        if (adjusted.model().initialize() != 0)
          return TestStatus(false).report(__func__);
        success *= adjusted.checkStateRows({{Internal::H, 1.0}}, "configured dam head after reinitialization");
        success *= adjusted.checkSteadyState();

        if (!adjusted.setPoint({.state = {{Internal::C, 0.75}, {Internal::RC, 0.2}}, .derivative = {{Internal::C, 0.0}}})
            || !adjusted.evaluateResidual())
          return TestStatus(false).report(__func__);
        success *= scalarMatches(
            adjusted.residual(Internal::C),
            0.0,
            "configured Gmax after reinitialization",
            kTol + 0.6 * std::exp(-0.2 * Math::MU<RealT>));

        // Initialization supports only a zero speed deviation; a moving
        // machine would need a multi-root gate search.
        success *= initializationRejectedAtomically(makeResidualData(),
                                                    0.4,
                                                    {{External::OMEGA, 0.03},
                                                     {External::PREF, 77.0},
                                                     {External::PAUX, 0.02}},
                                                    "nonzero initial speed deviation");

        // An invalid configuration is rejected before any state is written.
        Fixture<ScalarT> invalid(withParameters(makeResidualData(), {{Params::Rtemp, 0.0}}), __func__, kTol);
        configureHygov(invalid);
        attachInputs(invalid);
        if (invalid.prepare() || !poisonState(invalid, 0.4))
          return TestStatus(false).report(__func__);
        const auto invalid_before = invalid.snapshot();
        if (invalid.model().initialize() == 0)
        {
          std::cout << "Expected initialization rejection: invalid configuration\n";
          success = false;
        }
        success *= invalid.checkUnchanged(invalid_before);

        // Zero mechanical power lands on an in-range root and initializes at rest.
        Fixture<ScalarT> zero_power(makeData(), __func__, kTol);
        configureHygov(zero_power);
        if (!zero_power.initialize({{Internal::PMECH, 0.0}}))
          return TestStatus(false).report(__func__);
        success *= scalarMatches(zero_power.state(Internal::C), 0.1, "zero-power desired gate", kTol + 2.0 * rampError(0.1));
        success *= scalarMatches(zero_power.state(Internal::G), 0.1, "zero-power gate", kTol + 2.0 * rampError(0.1));
        success *= zero_power.checkSteadyState();

        // The smooth identity curve leaves a ln(2)/MU knee at each end, so
        // makeData()'s achievable component-base power range is
        // [knee - 0.1, 0.9 - knee].
        const RealT p_max = static_cast<RealT>(0.9) - knee;
        const RealT p_min = knee - static_cast<RealT>(0.1);

        Fixture<ScalarT> lower_edge(makeData(), __func__, kTol);
        configureHygov(lower_edge);
        if (!lower_edge.initialize({{Internal::PMECH, p_min - 0.5 * kTol}}))
          return TestStatus(false).report(__func__);
        success *= lower_edge.checkStateRows({{Internal::C, 0.0}, {Internal::G, 0.0}}, "half the tolerance below the achievable minimum");
        success *= scalarMatches(lower_edge.output(Internal::PMECH),
                                 p_min - 0.5 * kTol,
                                 "clipped pmech value");
        success *= lower_edge.checkSteadyState();

        Fixture<ScalarT> effective_edge(makeData(), __func__, kTol);
        configureHygov(effective_edge);
        if (!effective_edge.initialize({{Internal::PMECH, p_max + 0.5 * kTol}}))
          return TestStatus(false).report(__func__);
        success                         *= effective_edge.checkStateRows({{Internal::C, 1.0}, {Internal::G, 1.0}}, "half the tolerance beyond the achievable maximum");
        const RealT effective_edge_head  = effective_edge.state(Internal::H);
        if (!(effective_edge_head > 1.0))
        {
          std::cout << "effective head was not raised above configured Hdam\n";
          success = false;
        }
        success *= effective_edge.checkSteadyState();

        success *= initializationRejectedAtomically(
            makeData(),
            p_min - 2.0 * kTol,
            {{External::OMEGA, 0.0}, {External::PREF, 77.0}, {External::PAUX, 0.02}},
            "twice the tolerance below the achievable minimum");

        const RealT nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT infinity = std::numeric_limits<RealT>::infinity();

        // A non-finite input is rejected atomically, NaN included: the
        // exact-preservation check states what a tolerance comparison of a
        // NaN input never could.
        const std::array<RealT, 3> nonfinite_inputs{{nan, infinity, -infinity}};

        for (const RealT value : nonfinite_inputs)
        {
          success *= initializationRejectedAtomically(makeData(),
                                                      0.4,
                                                      {{External::OMEGA, value},
                                                       {External::PREF, 77.0},
                                                       {External::PAUX, 0.02}},
                                                      "non-finite speed input");
          success *= initializationRejectedAtomically(makeData(),
                                                      0.4,
                                                      {{External::OMEGA, 0.0},
                                                       {External::PREF, 77.0},
                                                       {External::PAUX, value}},
                                                      "non-finite auxiliary-power input");

          Fixture<ScalarT> pmech(makeData(), "non-finite pmech seed", kTol);
          configureHygov(pmech);
          attachInputs(pmech);
          if (!pmech.prepare()
              || !pmech.setPoint({.inputs = {{External::PREF, 77.0}, {External::PAUX, 0.02}},
                                  .state  = {{Internal::PMECH, value}}}))
            return TestStatus(false).report(__func__);
          const auto before  = pmech.snapshot();
          success           *= pmech.model().initialize() != 0;
          success           *= pmech.checkUnchanged(before);
        }

        return success.report(__func__);
      }

      /// Inversion may shift the gate near a curve corner; equilibrium must still be exact.
      TestOutcome initializationExactness()
      {
        TestStatus success = true;

        // At unit head, gain, and power base, these are the gate curve's own points.
        const std::array<std::pair<RealT, RealT>, 5> points{{
            {0.2, 0.15},
            {0.4, 0.42},
            {0.5, 0.54},
            {0.6, 0.66},
            {0.8, 0.85},
        }};
        for (const auto& [gate, power] : points)
        {
          Fixture<ScalarT> fixture(makeCurveData(), __func__, kTol);
          if (!initializeHygov(fixture, power))
            return TestStatus(false).report(__func__);
          // Inverting a curve with minimum slope 0.75 amplifies its value error.
          const RealT tolerance  = curveError(0.0) / 0.75;
          success               *= scalarMatches(fixture.state(Internal::C), gate, "desired gate", tolerance);
          success               *= scalarMatches(fixture.state(Internal::G), gate, "gate", tolerance);
          success               *= scalarMatches(fixture.output(Internal::PMECH), power, "preserved power");
          success               *= fixture.checkSteadyState();
        }
        return success.report(__func__);
      }

      /// Check all twelve equations against arithmetic and ideal limiter values.
      TestOutcome residualEquations()
      {
        Fixture<ScalarT> fixture(makeResidualData(), __func__, kTol);
        if (!initializeHygov(fixture, 0.4) || !fixture.setPoint(residualPoint()))
          return TestStatus(false).report(__func__);
        TestStatus success = fixture.checkResidualRows({
            {Internal::XN, -109.0 / 1400.0},
            {Internal::XF, -0.73},
            {Internal::G, 37.0 / 300.0},
            {Internal::Q, 3.0 / 260.0},
            {Internal::EF, 0.5863},
            {Internal::FC, -0.7405},
            {Internal::H, -0.0333},
            {Internal::PMECH, -0.01268},
        });

        success *= scalarMatches(fixture.residual(Internal::C), 0.06, "interior gate rate", kTol + 0.27 * std::exp(-0.43 * Math::MU<RealT>));
        success *= scalarMatches(fixture.residual(Internal::OMEGADB), 0.005, "speed outside the deadband", kTol + 0.04 * std::exp(-0.01 * Math::MU<RealT>));
        success *= scalarMatches(fixture.residual(Internal::RC), 0.03, "interior velocity", kTol + 2.0 * rampError(0.03));
        success *= scalarMatches(fixture.residual(Internal::PGV), -0.046, "gate-power interpolation", curveError(0.07));
        return success.report(__func__);
      }

      /// Speed deadband, desired-gate velocity limiting, and gate-position
      /// anti-windup.
      TestOutcome governorControl()
      {
        TestStatus success = true;
        const auto data    = makeResidualData();

        const RealT                       band          = 4.0 / Math::MU<RealT>;
        const auto                        deadband_data = withParameters(data, {{Params::db1, band}});
        const std::array<ResidualCase, 3> deadband_cases{{
            {"below the deadband",
             {.inputs = {{External::OMEGA, -3.0 * band}}, .state = {{Internal::OMEGADB, 0.0}}},
             {{Internal::OMEGADB, -3.0 * band}}},
            {"inside the deadband",
             {.inputs = {{External::OMEGA, band / 4.0}}, .state = {{Internal::OMEGADB, 0.0}}},
             {{Internal::OMEGADB, 0.0}}},
            {"above the deadband",
             {.inputs = {{External::OMEGA, 3.0 * band}}, .state = {{Internal::OMEGADB, 0.0}}},
             {{Internal::OMEGADB, 3.0 * band}}},
        }};
        // Sigmoid tails are bounded by exp(-MU * distance) on each side.
        for (const auto& test_case : deadband_cases)
        {
          const RealT      omega    = test_case.point.inputs.front().second;
          const RealT      distance = std::abs(std::abs(omega) - band);
          Fixture<ScalarT> fixture(deadband_data, test_case.label, kTol + 2.0 * std::abs(omega) * std::exp(-Math::MU<RealT> * distance));
          success *= initializeHygov(fixture, 0.4)
                     && fixture.checkResidualRows(test_case.point, test_case.expected);
        }

        const std::array<ResidualCase, 3> gate_velocity_cases{{
            {"gate velocity below the rate limit",
             {.state = {{Internal::FC, -0.6}, {Internal::RC, 0.0}}},
             {{Internal::RC, -0.15}}},
            {"gate velocity inside the rate limit",
             {.state = {{Internal::FC, 0.05}, {Internal::RC, 0.0}}},
             {{Internal::RC, 0.05}}},
            {"gate velocity above the rate limit",
             {.state = {{Internal::FC, 0.6}, {Internal::RC, 0.0}}},
             {{Internal::RC, 0.15}}},
        }};
        success *= runResidualCases(data, 0.4, gate_velocity_cases, kTol + 2.0 * rampError(0.1));

        const std::array<ResidualCase, 4> gate_antiwindup_cases{{
            {"Gmax blocks an outward desired-gate rate",
             {.state      = {{Internal::C, 1.2}, {Internal::RC, 0.2}},
              .derivative = {{Internal::C, 0.0}}},
             {{Internal::C, 0.0}}},
            {"Gmin blocks an outward desired-gate rate",
             {.state      = {{Internal::C, -0.2}, {Internal::RC, -0.2}},
              .derivative = {{Internal::C, 0.0}}},
             {{Internal::C, 0.0}}},
            {"Gmax admits a restoring desired-gate rate",
             {.state      = {{Internal::C, 1.2}, {Internal::RC, -0.2}},
              .derivative = {{Internal::C, 0.0}}},
             {{Internal::C, -0.2}}},
            {"Gmin admits a restoring desired-gate rate",
             {.state      = {{Internal::C, -0.2}, {Internal::RC, 0.2}},
              .derivative = {{Internal::C, 0.0}}},
             {{Internal::C, 0.2}}},
        }};
        success *= runResidualCases(data, 0.4, gate_antiwindup_cases, kTol + 0.6 * std::exp(-0.2 * Math::MU<RealT>));

        // The blocked gate retains F_yp = -1 and an explicit zero RC entry.
        const RealT                           margin = 40.0 / Math::MU<RealT>;
        Fixture<DependencyTracking::Variable> blocked(data, "blocked desired gate", kTol);
        success *= initializeHygov(blocked, 0.4)
                   && blocked.setPoint({.state      = {{Internal::C, 0.95 + margin}, {Internal::RC, margin}},
                                        .derivative = {{Internal::C, 0.0}}})
                   && blocked.checkJacobianRow(
                       Internal::C, {{Internal::C, -1.0}, {Internal::RC, 0.0}}, 1.0);

        return success.report(__func__);
      }

      /// Gate-power, water-column, damping, and curve-inversion behavior,
      /// including a flat segment.
      TestOutcome turbineDynamics()
      {
        TestStatus success = true;
        const auto data    = makeResidualData();

        const std::array<ResidualCase, 5> gate_power_cases{{
            {"gate-power curve segment 1",
             {.state = {{Internal::G, 0.1}, {Internal::PGV, 0.0}}},
             {{Internal::PGV, 0.075}}},
            {"gate-power curve segment 2",
             {.state = {{Internal::G, 0.3}, {Internal::PGV, 0.0}}},
             {{Internal::PGV, 0.285}}},
            {"gate-power curve segment 3",
             {.state = {{Internal::G, 0.5}, {Internal::PGV, 0.0}}},
             {{Internal::PGV, 0.54}}},
            {"gate-power curve segment 4",
             {.state = {{Internal::G, 0.7}, {Internal::PGV, 0.0}}},
             {{Internal::PGV, 0.755}}},
            {"gate-power curve segment 5",
             {.state = {{Internal::G, 0.9}, {Internal::PGV, 0.0}}},
             {{Internal::PGV, 0.925}}},
        }};
        success *= runResidualCases(data, 0.4, gate_power_cases, curveError(0.1));

        // A head away from the dam head drives the flow and head rows, and
        // turbine damping scales with speed deviation and gate.
        const std::array<ResidualCase, 2> turbine_cases{{
            {"water column",
             {.state      = {{Internal::Q, 0.61}, {Internal::H, 0.9}, {Internal::PGV, 0.55}},
              .derivative = {{Internal::Q, 0.05}}},
             {{Internal::Q, 47.0 / 260.0}, {Internal::H, -0.09985}}},
            {"turbine damping",
             {.inputs = {{External::OMEGA, 0.05}},
              .state  = {
                  {Internal::G, 0.6},
                  {Internal::Q, 0.7},
                  {Internal::H, 1.1},
                  {Internal::PMECH, 0.5},
              }},
             {{Internal::PMECH, -0.2678}}},
        }};
        success *= runResidualCases(data, 0.4, turbine_cases);

        Fixture<ScalarT> curve(makeCurveData(), __func__, kTol);
        if (!initializeHygov(curve, 0.54))
          return TestStatus(false).report(__func__);
        const RealT gate_error  = 2.0 * curveError(0.1) / 0.75;
        success                *= scalarMatches(curve.state(Internal::C), 0.5, "midpoint desired gate", gate_error);
        success                *= scalarMatches(curve.state(Internal::G), 0.5, "midpoint gate", gate_error);
        success                *= scalarMatches(curve.input(External::PREF), 0.03, "published pref", kTol + 0.06 * gate_error);
        success                *= scalarMatches(curve.output(Internal::PMECH), 0.54, "preserved power");
        success                *= curve.checkSteadyState();

        // A flat source-curve segment must initialize to a gate on that segment.
        // makeData() uses equal power bases, At = Hdam = 1, and Qnl = 0.1,
        // so a 0.5 plateau maps to pmech = 0.4 without encoding Math::MU.
        const RealT      flat_gate_minimum = static_cast<RealT>(0.4);
        const RealT      flat_gate_maximum = static_cast<RealT>(0.6);
        const RealT      plateau_power     = static_cast<RealT>(0.5);
        const RealT      plateau_pmech     = static_cast<RealT>(0.4);
        Fixture<ScalarT> flat(withParameters(makeData(), {{Params::Pgv2, plateau_power}, {Params::Pgv3, plateau_power}}), __func__, kTol);
        configureHygov(flat);
        if (!flat.initialize({{Internal::PMECH, plateau_pmech}}))
          return TestStatus(false).report(__func__);
        const RealT flat_gate =
            flat.state(Internal::G);
        if (flat_gate < flat_gate_minimum || flat_gate > flat_gate_maximum)
        {
          std::cout << "flat-segment plateau gate "
                    << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                    << flat_gate << " is outside [" << flat_gate_minimum
                    << ", " << flat_gate_maximum << "]\n";
          success = false;
        }
        success *= flat.checkSteadyState();

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Every Enzyme CSR row must match dependency tracking at gates inside
      /// each curve segment and at each breakpoint, and both paths must
      /// carry the PGV row's gate dependence.
      TestOutcome jacobian()
      {
        TestStatus success = true;
        const auto data    = makeResidualData();
        for (const RealT gate : {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9})
        {
          const std::string label = "HYGOV gate " + std::to_string(gate);
          const auto        setup = [&](auto& fixture)
          {
            return initializeHygov(fixture, 0.4)
                   && fixture.setPoint(residualPoint())
                   && fixture.setState({{Internal::G, gate}});
          };
          success *= Fixture<ScalarT>::checkJacobian(
              data, setup, {0.0, 1.0, 2.5}, label.c_str(), kTol, {{Internal::PGV, Internal::G}});
        }
        return success.report(__func__);
      }
#endif

    private:
      using Params   = PhasorDynamics::Governor::HygovParameters;
      using Internal = PhasorDynamics::Governor::HygovInternalVariables;
      using External = PhasorDynamics::Governor::HygovExternalVariables;
      using Mon      = PhasorDynamics::Governor::HygovMonitorableVariables;
      using Data     = PhasorDynamics::Governor::HygovData<RealT, IdxT>;
      using HygovT   = PhasorDynamics::Governor::Hygov<ScalarT, IdxT>;

      template <typename T>
      using Fixture = ComponentTestFixture<PhasorDynamics::Governor::Hygov, T, IdxT>;
      using Point   = typename Fixture<ScalarT>::Point;
      using Values  = typename Fixture<ScalarT>::Values;
      using Inputs  = typename Fixture<ScalarT>::Inputs;

      struct ResidualCase
      {
        const char* label;
        Point       point;
        Values      expected;
      };

      // Softplus differs from the ideal ramp by at most exp(-MU * distance) / MU.
      static RealT rampError(RealT distance)
      {
        return std::exp(-Math::MU<RealT> * distance) / Math::MU<RealT>;
      }

      // The nonidentity curve's absolute slope changes sum to 2.7.
      static RealT curveError(RealT distance)
      {
        return kTol + 2.7 * rampError(distance);
      }

      static Data withParameters(Data                                            data,
                                 std::initializer_list<std::pair<Params, RealT>> overrides)
      {
        for (const auto& [parameter, value] : overrides)
        {
          data.parameters[parameter] = value;
        }
        return data;
      }

      template <typename T>
      void configureHygov(Fixture<T>& fixture, RealT system_va_base = 100.0e6) const
      {
        fixture.model().setSystemBase(60.0, system_va_base);
        fixture.template assignOutput<Internal::PMECH>();
      }

      template <typename T>
      void attachInputs(Fixture<T>& fixture) const
      {
        fixture.template attachInput<External::OMEGA>(0.0);
        fixture.template attachInput<External::PREF>(0.0);
        fixture.template attachInput<External::PAUX>(0.0);
      }

      template <typename T>
      bool initializeHygov(Fixture<T>& fixture, RealT pmech) const
      {
        configureHygov(fixture);
        attachInputs(fixture);
        return fixture.initialize({{Internal::PMECH, pmech}});
      }

      struct RestoreVerbosity
      {
        const Log::Verbosity previous = Log::verbosity();

        ~RestoreVerbosity()
        {
          Log::setVerbosity(previous);
        }
      };

      Data makeMinimalData() const
      {
        Data data;
        data.device_class              = "Hygov";
        data.disambiguation_string     = "hygov_test";
        data.parameters[Params::Trate] = 100.0;
        data.monitored_variables.insert(Mon::pmech);
        data.monitored_variables.insert(Mon::filter);
        data.monitored_variables.insert(Mon::desiredgate);
        data.monitored_variables.insert(Mon::gate);
        data.monitored_variables.insert(Mon::flow);
        data.monitored_variables.insert(Mon::head);
        return data;
      }

      Data makeExplicitDefaultData() const
      {
        // These are the documented defaults. The all-zero source curve
        // selects the identity curve, spelled out here point by point.
        return withParameters(makeMinimalData(),
                              {{Params::Rperm, 0.04},
                               {Params::Rtemp, 0.3},
                               {Params::Tr, 5.0},
                               {Params::Tf, 0.05},
                               {Params::Tg, 0.5},
                               {Params::Velm, 0.2},
                               {Params::Gmax, 1.0},
                               {Params::Gmin, 0.0},
                               {Params::Tw, 1.0},
                               {Params::At, 1.2},
                               {Params::Dturb, 0.5},
                               {Params::Qnl, 0.05},
                               {Params::Tn, 0.0},
                               {Params::Tnp, 0.0},
                               {Params::db1, 0.0},
                               {Params::db2, 0.0},
                               {Params::Hdam, 1.0},
                               {Params::Gv0, 0.0},
                               {Params::Gv1, 0.2},
                               {Params::Gv2, 0.4},
                               {Params::Gv3, 0.6},
                               {Params::Gv4, 0.8},
                               {Params::Gv5, 1.0},
                               {Params::Pgv0, 0.0},
                               {Params::Pgv1, 0.2},
                               {Params::Pgv2, 0.4},
                               {Params::Pgv3, 0.6},
                               {Params::Pgv4, 0.8},
                               {Params::Pgv5, 1.0}});
      }

      Data makeData() const
      {
        // The documented typical values with the floored time constants
        // raised above the floor, so routine fixtures log no warnings.
        return withParameters(makeMinimalData(),
                              {{Params::Trate, 100.0},
                               {Params::Rperm, 0.05},
                               {Params::Rtemp, 0.4},
                               {Params::Tr, 5.0},
                               {Params::Tf, 0.2},
                               {Params::Tg, 0.5},
                               {Params::Velm, 0.5},
                               {Params::Gmax, 1.0},
                               {Params::Gmin, 0.0},
                               {Params::Tw, 1.0},
                               {Params::At, 1.0},
                               {Params::Dturb, 0.0},
                               {Params::Qnl, 0.1},
                               {Params::Tn, 0.0},
                               {Params::Tnp, 1.0},
                               {Params::db1, 0.0},
                               {Params::db2, 0.0},
                               {Params::Hdam, 1.0},
                               {Params::Gv0, 0.0},
                               {Params::Gv1, 0.2},
                               {Params::Gv2, 0.4},
                               {Params::Gv3, 0.6},
                               {Params::Gv4, 0.8},
                               {Params::Gv5, 1.0},
                               {Params::Pgv0, 0.0},
                               {Params::Pgv1, 0.2},
                               {Params::Pgv2, 0.4},
                               {Params::Pgv3, 0.6},
                               {Params::Pgv4, 0.8},
                               {Params::Pgv5, 1.0}});
      }

      Data makeResidualData() const
      {
        return withParameters(makeData(),
                              {{Params::Trate, 50.0},
                               {Params::Rperm, 0.06},
                               {Params::Rtemp, 0.4},
                               {Params::Tr, 4.0},
                               {Params::Tf, 0.2},
                               {Params::Tg, 0.6},
                               {Params::Velm, 0.15},
                               {Params::Gmax, 0.95},
                               {Params::Gmin, 0.05},
                               {Params::Tw, 1.3},
                               {Params::At, 1.1},
                               {Params::Dturb, 0.6},
                               {Params::Qnl, 0.08},
                               {Params::Tn, 0.7},
                               {Params::Tnp, 1.4},
                               {Params::db1, 0.01},
                               {Params::Hdam, 1.2},
                               {Params::Pgv1, 0.15},
                               {Params::Pgv2, 0.42},
                               {Params::Pgv3, 0.66},
                               {Params::Pgv4, 0.85}});
      }

      Data makeCurveData() const
      {
        return withParameters(makeResidualData(),
                              {{Params::Trate, 100.0}, {Params::At, 1.0}, {Params::Qnl, 0.0}, {Params::Hdam, 1.0}});
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      Point residualPoint() const
      {
        return {
            .inputs = {{External::OMEGA, 0.02}, {External::PREF, 0.31}, {External::PAUX, 0.07}},
            .state  = {
                {Internal::XN, 0.11},
                {Internal::XF, 0.23},
                {Internal::C, 0.52},
                {Internal::G, 0.47},
                {Internal::Q, 0.61},
                {Internal::OMEGADB, 0.015},
                {Internal::EF, 0.08},
                {Internal::FC, 0.12},
                {Internal::RC, 0.09},
                {Internal::PGV, 0.55},
                {Internal::H, 1.12},
                {Internal::PMECH, 0.33},
            },
            .derivative = {
                {Internal::XN, 0.01},
                {Internal::XF, -0.02},
                {Internal::C, 0.03},
                {Internal::G, -0.04},
                {Internal::Q, 0.05},
            }};
      }

      /// Omitting every optional parameter must give exactly the model built
      /// from the defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), __func__, kTol);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), __func__, kTol);
        configureHygov(implicit_defaults, 200.0e6);
        configureHygov(explicit_defaults, 200.0e6);
        attachInputs(implicit_defaults);
        attachInputs(explicit_defaults);
        if (!implicit_defaults.initialize({{Internal::PMECH, 0.3}})
            || !explicit_defaults.initialize({{Internal::PMECH, 0.3}})
            || !explicit_defaults.evaluateResidual())
          return false;
        const auto reference  = explicit_defaults.snapshot();
        bool       success    = implicit_defaults.checkStateRows(reference.state);
        success              &= implicit_defaults.checkDerivativeRows(reference.derivative);
        success              &= implicit_defaults.checkResiduals(explicit_defaults.residuals());
        if (!explicit_defaults.setPoint(residualPoint()) || !explicit_defaults.evaluateResidual())
          return false;
        success &= implicit_defaults.checkResiduals(residualPoint(), explicit_defaults.residuals());
        return success;
      }

      template <External variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeData(), __func__, kTol);
        configureHygov(fixture);
        fixture.model().getSignals().template attachSignalNode<variable>(&unlinked_node);
        return fixture.model().verify() > 0;
      }

      /// Fill the state and derivative with a recognizable ramp, then restore
      /// the aliased pmech entry, so any write by a rejected initialization
      /// is visible.
      bool poisonState(Fixture<ScalarT>& fixture, RealT pmech) const
      {
        auto* y  = fixture.model().y().getData();
        auto* yp = fixture.model().yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.model().y().getSize()); ++i)
        {
          y[i]  = 0.125 + 0.01 * static_cast<RealT>(i);
          yp[i] = -0.25 - 0.01 * static_cast<RealT>(i);
        }
        fixture.model().y().setDataUpdated();
        fixture.model().yp().setDataUpdated();
        return fixture.setState({{Internal::PMECH, pmech}});
      }

      /// Initialization must fail and leave the poisoned state, the seeded
      /// pmech value, and every supplied input untouched.
      bool initializationRejectedAtomically(const Data& data, RealT pmech, const Inputs& inputs, const char* label) const
      {
        Fixture<ScalarT> fixture(data, label, kTol);
        configureHygov(fixture);
        attachInputs(fixture);
        if (!fixture.prepare() || !fixture.setPoint({.inputs = inputs})
            || !poisonState(fixture, pmech))
          return false;
        const auto before   = fixture.snapshot();
        const bool rejected = fixture.model().initialize() != 0;
        if (!rejected)
          std::cout << "Expected initialization rejection: " << label << '\n';
        const bool unchanged = fixture.checkUnchanged(before);
        return rejected && unchanged;
      }

      template <size_t size>
      bool runResidualCases(const Data& data, RealT pmech, const std::array<ResidualCase, size>& cases, RealT tolerance = kTol) const
      {
        bool success = true;
        for (const auto& test_case : cases)
        {
          Fixture<ScalarT> fixture(data, test_case.label, tolerance);
          success &= initializeHygov(fixture, pmech)
                     && fixture.checkResidualRows(
                         test_case.point,
                         test_case.expected);
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
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << "\n";
        return false;
      }
    };
  } // namespace Testing
} // namespace GridKit
