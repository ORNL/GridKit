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
#include <GridKit/Utilities/MapFromCsr.hpp>

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

        Fixture<ScalarT> configured(makeData());
        success *= (configured.hygov.size() == static_cast<IdxT>(Internal::MAXIMUM));
        success *= (configured.hygov.getMonitor() != nullptr);
        success *= (configured.hygov.verify() == 0);

        noteExpectedLogs("Testing HYGOV defaults and invalid configurations. "
                         "Logged errors, time-constant warnings, and an unsupported "
                         "backlash warning are expected.");

        Fixture<ScalarT> minimal(makeMinimalData());
        success *= (minimal.hygov.verify() == 0);
        success *= defaultsMatchDocumentedValues();

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
            Fixture<ScalarT> invalid_fixture(makeData(), {{parameter, value}});
            success *= (invalid_fixture.hygov.verify() > 0);
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
          Fixture<ScalarT> invalid_fixture(makeData(), {{parameter, value}});
          success *= (invalid_fixture.hygov.verify() > 0);
        }

        // A curve with no rise cannot yield a unique gate.
        Fixture<ScalarT> flat_curve(makeData(),
                                    {{Params::Pgv1, 0.0},
                                     {Params::Pgv2, 0.0},
                                     {Params::Pgv3, 0.0},
                                     {Params::Pgv4, 0.0},
                                     {Params::Pgv5, 0.0}});
        success *= (flat_curve.hygov.verify() > 0);

        // A curve that rises only outside the configured response limits is
        // valid because initialization may expand those limits.
        Fixture<ScalarT> flat_configured_range(
            makeData(),
            {{Params::Gmin, 0.0},
             {Params::Gmax, 0.2},
             {Params::Pgv0, 0.5},
             {Params::Pgv1, 0.5},
             {Params::Pgv2, 0.5},
             {Params::Pgv3, 0.5},
             {Params::Pgv4, 0.5},
             {Params::Pgv5, 1.0}});
        success *= (flat_configured_range.hygov.verify() == 0);

        // A requested backlash is accepted, warns, and remains inactive.
        Fixture<ScalarT> backlash(makeData(), {{Params::db2, 0.5}});
        success *= (backlash.hygov.verify() == 0);

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                   = makeData();
        integer_real.parameters[Params::Tw] = static_cast<IdxT>(2);
        Fixture<ScalarT> integer_model(integer_real);
        success *= (integer_model.hygov.verify() == 0);

        auto bad_numeric_type                      = makeData();
        bad_numeric_type.parameters[Params::Trate] = true;
        Fixture<ScalarT> bad_numeric_model(bad_numeric_type);
        success *= (bad_numeric_model.hygov.verify() > 0);

        Fixture<ScalarT> overflowing_component_base(
            makeData(),
            {{Params::Trate, std::numeric_limits<RealT>::max()}});
        success *= (overflowing_component_base.hygov.verify() > 0);

        Fixture<ScalarT> overflowing_base_ratio(
            makeData(),
            {{Params::Trate, std::numeric_limits<RealT>::min()}});
        success *= (overflowing_base_ratio.hygov.verify() > 0);

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
          Fixture<ScalarT> invalid_base(makeData(), {}, system_base);
          success *= (invalid_base.hygov.verify() > 0);
        }

        success *= unlinkedSignalRejected<External::OMEGA>();
        success *= unlinkedSignalRejected<External::PREF>();
        success *= unlinkedSignalRejected<External::PAUX>();

        // All five zero time constants use the documented numerical floor and
        // still admit a consistent steady-state initialization.
        Fixture<ScalarT> floors(makeData(),
                                {{Params::Tr, 0.0},
                                 {Params::Tf, 0.0},
                                 {Params::Tg, 0.0},
                                 {Params::Tw, 0.0},
                                 {Params::Tnp, 0.0}});
        success *= floors.initialize(0.4);
        success *= (floors.evaluate() == 0);
        success *= allResidualsZero(floors.hygov);

        return success.report(__func__);
      }

      /// A nonidentity power-base initialization with every port attached.
      /// The machine-provided pmech value must remain unchanged while HYGOV
      /// initializes and publishes its resolved load reference.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeData(), {{Params::Trate, 50.0}});
        fixture.attachAllInputs();
        fixture.input(External::PAUX)  = 0.02;
        fixture.input(External::PREF)  = 99.0; // stale value the publication must replace
        success                       *= fixture.initialize(0.4);
        success                       *= (fixture.hygov.tagDifferentiable() == 0);
        success                       *= (fixture.evaluate() == 0);

        const auto* y  = fixture.hygov.y().getData();
        success       *= scalarMatches(y[static_cast<size_t>(Internal::XF)], 0.0, "XF at rest");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::C)],
                                 0.9000000000001573,
                                 "C on component base");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::G)],
                                 0.9000000000001573,
                                 "G on component base");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::Q)], 0.9, "Q on component base");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::PGV)],
                                 0.9,
                                 "PGV on component base");
        success       *= scalarMatches(y[static_cast<size_t>(Internal::H)], 1.0, "H at the dam head");
        success       *= scalarMatches(fixture.pmech(), 0.4, "preserved pmech value");

        success *= scalarMatches(fixture.input(External::OMEGA), 0.0, "preserved omega input");
        success *= scalarMatches(fixture.input(External::PREF), 0.0025, "published pref");
        success *= scalarMatches(fixture.input(External::PAUX), 0.02, "preserved paux input");

        // Verify the six documented outputs through the public monitor controller.
        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.hygov.getMonitor());
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
          success *= scalarMatches(monitored[3], 0.9000000000001573, "monitored desiredgate");
          success *= scalarMatches(monitored[4], 0.9000000000001573, "monitored gate");
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
        for (size_t i = 0; i < static_cast<size_t>(fixture.hygov.size()); ++i)
        {
          const bool differential = i <= static_cast<size_t>(Internal::Q);
          if (fixture.hygov.tag()[i] != differential)
          {
            std::cout << "HYGOV differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }

        success *= allResidualsZero(fixture.hygov);

        // A system-base reference step lands on the governor error scaled by
        // the base ratio.
        fixture.input(External::PREF)  = 0.1025; // the published 0.0025 plus a 0.1 step
        success                       *= (fixture.evaluate() == 0);
        success                       *= residualsMatch(fixture.hygov,
                                                        {{Internal::EF, 0.2}},
                                  "reference step on the component base");

        // Unattached ports fall back to the references latched by
        // initialize(), so the same steady state holds without a controller.
        Fixture<ScalarT> fallback(makeData(), {{Params::Trate, 50.0}});
        success *= fallback.initialize(0.4);
        success *= (fallback.evaluate() == 0);
        success *= allResidualsZero(fallback.hygov);

        return success.report(__func__);
      }

      /// Mechanical-power, gate-limit, speed-deviation, and finite-input
      /// initialization domains. Response limits and high-power dam head
      /// are adjusted when needed; every rejected initialization is atomic.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing HYGOV initialization boundaries. "
                         "Logged errors, response-limit warnings, and dam-head warnings "
                         "are expected.");

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

        // 4.5 MW on the system base is 2.5 pu on a 1.8 MW turbine base.
        Fixture<ScalarT> effective_fixture(
            makeData(),
            {{Params::Trate, 1.8},
             {Params::At, 1.25},
             {Params::Qnl, 0.07},
             {Params::Gmax, 0.5}});
        effective_fixture.attachAllInputs();
        success *= effective_fixture.initialize(0.045);
        success *= stateMatches(
            effective_fixture.hygov,
            {{Internal::C, 1.0},
             {Internal::G, 1.0},
             {Internal::Q, 1.2812656647316965},
             {Internal::PGV, 0.9971118867476669},
             {Internal::H, 1.6511654364800423}},
            "effective dam head");
        success *= scalarMatches(effective_fixture.pmech(), 0.045, "preserved pmech value");
        success *= scalarMatches(effective_fixture.input(External::PREF),
                                 0.0009,
                                 "published pref");
        success *= (effective_fixture.evaluate() == 0);
        success *= allResidualsZero(effective_fixture.hygov);

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
          Fixture<ScalarT> fixture(makeResidualData(),
                                   {{test_case.limit_parameter, test_case.limit}});
          success          *= fixture.initialize(0.4);
          const RealT gate  = static_cast<RealT>(
              fixture.hygov.y().getData()[static_cast<size_t>(Internal::C)]);
          const bool gate_is_outside = test_case.rate > 0.0
                                           ? gate > test_case.limit
                                           : gate < test_case.limit;
          if (!gate_is_outside)
          {
            std::cout << test_case.label << " did not initialize outside the configured limit\n";
            success = false;
          }
          success *= stateMatches(fixture.hygov,
                                  {{Internal::G, gate}, {Internal::H, 1.2}},
                                  test_case.label);
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.hygov);

          // The effective response bound admits an outward rate between the
          // configured limit and initialized gate.
          setState(fixture.hygov,
                   {{Internal::C, 0.5 * (test_case.limit + gate)},
                    {Internal::RC, test_case.rate}});
          setDerivative(fixture.hygov, {{Internal::C, 0.0}});
          success                   *= (fixture.evaluate() == 0);
          const RealT response_rate  = static_cast<RealT>(
              fixture.hygov.getResidual().getData()[static_cast<size_t>(Internal::C)]);
          const bool rate_is_admitted = test_case.rate > 0.0
                                            ? response_rate > 0.9 * test_case.rate
                                            : response_rate < 0.9 * test_case.rate;
          if (!rate_is_admitted)
          {
            std::cout << test_case.label << " did not admit the outward desired-gate rate\n";
            success = false;
          }
        }

        // A failed retry preserves the effective head and response bounds from
        // the prior success.
        const auto effective_y                    = copyVector(effective_fixture.hygov.y());
        const auto effective_yp                   = copyVector(effective_fixture.hygov.yp());
        effective_fixture.input(External::OMEGA)  = 0.03;
        success                                  *= (effective_fixture.hygov.initialize() != 0);
        success                                  *= vectorUnchanged(effective_fixture.hygov.y(), effective_y, "state");
        success                                  *= vectorUnchanged(effective_fixture.hygov.yp(), effective_yp, "derivative");
        success                                  *= scalarMatches(effective_fixture.input(External::PREF),
                                 0.0009,
                                 "preserved pref");
        effective_fixture.input(External::OMEGA)  = 0.0;
        success                                  *= (effective_fixture.evaluate() == 0);
        success                                  *= allResidualsZero(effective_fixture.hygov);

        setState(effective_fixture.hygov,
                 {{Internal::C, 0.75}, {Internal::RC, 0.2}});
        setDerivative(effective_fixture.hygov, {{Internal::C, 0.0}});
        success                    *= (effective_fixture.evaluate() == 0);
        const RealT preserved_rate  = static_cast<RealT>(
            effective_fixture.hygov.getResidual().getData()[static_cast<size_t>(Internal::C)]);
        if (!(preserved_rate > 0.19))
        {
          std::cout << "failed initialization did not preserve effective Gmax\n";
          success = false;
        }

        // A later feasible initialization starts again from configured limits
        // and Hdam.
        effective_fixture.setPmech(0.009);
        success *= (effective_fixture.hygov.initialize() == 0);
        success *= stateMatches(effective_fixture.hygov,
                                {{Internal::H, 1.0}},
                                "configured dam head after reinitialization");
        success *= (effective_fixture.evaluate() == 0);
        success *= allResidualsZero(effective_fixture.hygov);

        setState(effective_fixture.hygov,
                 {{Internal::C, 0.75}, {Internal::RC, 0.2}});
        setDerivative(effective_fixture.hygov, {{Internal::C, 0.0}});
        success *= (effective_fixture.evaluate() == 0);
        success *= scalarMatches(
            static_cast<RealT>(effective_fixture.hygov.getResidual().getData()[static_cast<size_t>(Internal::C)]),
            0.0,
            "configured Gmax after reinitialization");

        // Initialization supports only a zero speed deviation; a moving
        // machine would need a multi-root gate search.
        success *= initializationRejectedAtomically(makeResidualData(),
                                                    0.4,
                                                    {{External::OMEGA, 0.03},
                                                     {External::PREF, 77.0},
                                                     {External::PAUX, 0.02}},
                                                    "nonzero initial speed deviation");

        // An invalid configuration is rejected before any state is written.
        Fixture<ScalarT> invalid_fixture(makeResidualData(), {{Params::Rtemp, 0.0}});
        invalid_fixture.attachAllInputs();
        success *= (invalid_fixture.hygov.allocate() == 0);
        poisonState(invalid_fixture, 0.4);
        const auto invalid_y  = copyVector(invalid_fixture.hygov.y());
        const auto invalid_yp = copyVector(invalid_fixture.hygov.yp());
        if (invalid_fixture.hygov.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: invalid configuration\n";
          success = false;
        }
        success *= vectorUnchanged(invalid_fixture.hygov.y(), invalid_y, "state");
        success *= vectorUnchanged(invalid_fixture.hygov.yp(), invalid_yp, "derivative");

        // Zero mechanical power lands on an in-range root and initializes at rest.
        Fixture<ScalarT> zero_power_fixture(makeData());
        success *= zero_power_fixture.initialize(0.0);
        success *= stateMatches(
            zero_power_fixture.hygov,
            {{Internal::C, 0.09999999999984271}, {Internal::G, 0.09999999999984271}},
            "zero mechanical power");
        success *= (zero_power_fixture.evaluate() == 0);
        success *= allResidualsZero(zero_power_fixture.hygov);

        // The smooth identity curve leaves a ln(2)/MU knee at each end, so
        // makeData()'s achievable component-base power range is
        // [knee - 0.1, 0.9 - knee].
        const RealT knee  = std::log(static_cast<RealT>(2.0)) / Math::MU<RealT>;
        const RealT p_max = static_cast<RealT>(0.9) - knee;
        const RealT p_min = knee - static_cast<RealT>(0.1);

        Fixture<ScalarT> lower_edge(makeData());
        success *= lower_edge.initialize(p_min - 0.5 * kTol);
        success *= stateMatches(lower_edge.hygov,
                                {{Internal::C, 0.0}, {Internal::G, 0.0}},
                                "half the tolerance below the achievable minimum");
        success *= scalarMatches(lower_edge.pmech(),
                                 p_min - 0.5 * kTol,
                                 "clipped pmech value");
        success *= (lower_edge.evaluate() == 0);
        success *= allResidualsZero(lower_edge.hygov);

        Fixture<ScalarT> effective_edge(makeData());
        success                         *= effective_edge.initialize(p_max + 0.5 * kTol);
        success                         *= stateMatches(effective_edge.hygov,
                                                        {{Internal::C, 1.0}, {Internal::G, 1.0}},
                                "half the tolerance beyond the achievable maximum");
        const RealT effective_edge_head  = static_cast<RealT>(
            effective_edge.hygov.y().getData()[static_cast<size_t>(Internal::H)]);
        if (!(effective_edge_head > 1.0))
        {
          std::cout << "effective head was not raised above configured Hdam\n";
          success = false;
        }
        success *= (effective_edge.evaluate() == 0);
        success *= allResidualsZero(effective_edge.hygov);

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

          // A non-finite seed lands in the aliased pmech state itself, so the
          // poisoned-state comparison cannot express its preservation. The
          // inputs still must survive untouched.
          Fixture<ScalarT> pmech_fixture(makeData());
          pmech_fixture.attachAllInputs();
          pmech_fixture.input(External::PREF)  = 77.0;
          pmech_fixture.input(External::PAUX)  = 0.02;
          success                             *= pmech_fixture.prepare(value);
          success                             *= (pmech_fixture.hygov.initialize() != 0);
          success                             *= scalarPreserved(
              static_cast<RealT>(pmech_fixture.input(External::PREF)),
              77.0,
              "external input",
              static_cast<size_t>(External::PREF));
          success *= scalarPreserved(static_cast<RealT>(pmech_fixture.input(External::PAUX)),
                                     0.02,
                                     "external input",
                                     static_cast<size_t>(External::PAUX));
        }

        return success.report(__func__);
      }

      /// Initialization solves the smooth gate curve the residual evaluates,
      /// so every steady residual rests at machine rounding even where the
      /// smoothing bends the curve away from its piecewise-linear points.
      TestOutcome initializationExactness()
      {
        TestStatus success = true;

        // Values landing mid-segment and within the smoothing knee of every
        // interior curve breakpoint, where a piecewise-linear inversion
        // misses the implemented curve by up to O(1e-3). The gate literal
        // proves where each value lands.
        struct ExactnessCase
        {
          const char* label;
          RealT       pmech;
          RealT       gate;
        };

        const std::array<ExactnessCase, 5> exactness_cases{{
            {"gate inside the Gv1 knee", 0.0556, 0.1982318164100278},
            {"gate inside the Gv2 knee", 0.2509, 0.4003865335541374},
            {"gate mid-segment", 0.4, 0.5719050089028755},
            {"gate inside the Gv3 knee", 0.4244, 0.6007061471851347},
            {"gate inside the Gv4 knee", 0.5617, 0.8006094230811988},
        }};

        for (const auto& seed : exactness_cases)
        {
          Fixture<ScalarT> fixture(makeResidualData());
          success *= fixture.initialize(seed.pmech);
          success *= stateMatches(fixture.hygov, {{Internal::G, seed.gate}}, seed.label);
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.hygov);
        }

        return success.report(__func__);
      }

      /// A fixed numerical answer key for all 12 HYGOV residual rows. The
      /// expected values are literals, not a second implementation of HYGOV.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.hygov);
        success *= (fixture.evaluate() == 0);

        const std::array<InternalRow, static_cast<size_t>(Internal::MAXIMUM)> expected{{
            {Internal::XN, -0.07785714285714286},
            {Internal::XF, -0.7300000000000001},
            {Internal::C, 0.06},
            {Internal::G, 0.1233333333333334},
            {Internal::Q, 0.011538461538461414},
            {Internal::OMEGADB, 0.0033514666467982894},
            {Internal::EF, 0.5863},
            {Internal::FC, -0.7405000000000002},
            {Internal::RC, 0.029996890386450745},
            {Internal::PGV, -0.04600000003160343},
            {Internal::H, -0.033299999999999885},
            {Internal::PMECH, -0.012679999999999934},
        }};

        success *= (static_cast<size_t>(fixture.hygov.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.hygov, expected, "answer key");

        return success.report(__func__);
      }

      /// Speed deadband, desired-gate velocity limiting, and gate-position
      /// anti-windup.
      TestOutcome governorControl()
      {
        TestStatus success = true;
        const auto data    = makeResidualData();

        // Exercise both sides and the interior of the type-1 +/-0.01 deadband.
        const std::array<ResidualCase, 3> deadband_cases{{
            {"speed deadband below the band",
             {{External::OMEGA, -0.05}},
             {{Internal::OMEGADB, 0.0}},
             {},
             {{Internal::OMEGADB, -0.049996641662021946}}},
            {"speed deadband inside the band",
             {{External::OMEGA, 0.004}},
             {{Internal::OMEGADB, 0.0}},
             {},
             {{Internal::OMEGADB, 0.0009004582873718001}}},
            {"speed deadband above the band",
             {{External::OMEGA, 0.05}},
             {{Internal::OMEGADB, 0.0}},
             {},
             {{Internal::OMEGADB, 0.049996641662021946}}},
        }};
        success *= runResidualCases(data, 0.4, deadband_cases);

        const std::array<ResidualCase, 3> gate_velocity_cases{{
            {"gate velocity below the rate limit",
             {},
             {{Internal::FC, -0.6}, {Internal::RC, 0.0}},
             {},
             {{Internal::RC, -0.15}}},
            {"gate velocity inside the rate limit",
             {},
             {{Internal::FC, 0.05}, {Internal::RC, 0.0}},
             {},
             {{Internal::RC, 0.04999999999984272}}},
            {"gate velocity above the rate limit",
             {},
             {{Internal::FC, 0.6}, {Internal::RC, 0.0}},
             {},
             {{Internal::RC, 0.15000000000000002}}},
        }};
        success *= runResidualCases(data, 0.4, gate_velocity_cases);

        const std::array<ResidualCase, 4> gate_antiwindup_cases{{
            {"Gmax blocks an outward desired-gate rate",
             {},
             {{Internal::C, 1.2}, {Internal::RC, 0.2}},
             {{Internal::C, 0.0}},
             {{Internal::C, 0.0}}},
            {"Gmin blocks an outward desired-gate rate",
             {},
             {{Internal::C, -0.2}, {Internal::RC, -0.2}},
             {{Internal::C, 0.0}},
             {{Internal::C, 0.0}}},
            {"Gmax admits a restoring desired-gate rate",
             {},
             {{Internal::C, 1.2}, {Internal::RC, -0.2}},
             {{Internal::C, 0.0}},
             {{Internal::C, -0.2}}},
            {"Gmin admits a restoring desired-gate rate",
             {},
             {{Internal::C, -0.2}, {Internal::RC, 0.2}},
             {{Internal::C, 0.0}},
             {{Internal::C, 0.2}}},
        }};
        success *= runResidualCases(data, 0.4, gate_antiwindup_cases);

        // At alpha = 1, a blocked desired-gate row has derivative coefficient
        // -1 and no RC dependence, independently of either Jacobian backend.
        {
          using DepVar = DependencyTracking::Variable;

          Fixture<DepVar> blocked(data);
          blocked.attachAllInputs();
          success *= blocked.initialize(0.4);
          setState(blocked.hygov, {{Internal::C, 1.2}, {Internal::RC, 0.2}});
          setDerivative(blocked.hygov, {{Internal::C, 0.0}});
          numberVariables(blocked);
          success *= (blocked.evaluate() == 0);

          const auto& dependencies =
              blocked.hygov.getResidual().getData()[static_cast<size_t>(Internal::C)].getDependencies();
          const DepVar::DependencyMap expected{{
              {static_cast<size_t>(Internal::C), -1.0},
              {static_cast<size_t>(Internal::RC), 0.0},
          }};
          success *= isEqual(dependencies, expected, kTol);
        }

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
             {},
             {{Internal::G, 0.1}, {Internal::PGV, 0.0}},
             {},
             {{Internal::PGV, 0.07500000000021236}}},
            {"gate-power curve segment 2",
             {},
             {{Internal::G, 0.3}, {Internal::PGV, 0.0}},
             {},
             {{Internal::PGV, 0.28500000000007075}}},
            {"gate-power curve segment 3",
             {},
             {{Internal::G, 0.5}, {Internal::PGV, 0.0}},
             {},
             {{Internal::PGV, 0.5399999999999371}}},
            {"gate-power curve segment 4",
             {},
             {{Internal::G, 0.7}, {Internal::PGV, 0.0}},
             {},
             {{Internal::PGV, 0.7549999999999292}}},
            {"gate-power curve segment 5",
             {},
             {{Internal::G, 0.9}, {Internal::PGV, 0.0}},
             {},
             {{Internal::PGV, 0.9249999999998506}}},
        }};
        success *= runResidualCases(data, 0.4, gate_power_cases);

        // A head away from the dam head drives the flow and head rows, and
        // turbine damping scales with speed deviation and gate.
        const std::array<ResidualCase, 2> turbine_cases{{
            {"water column",
             {},
             {{Internal::Q, 0.61}, {Internal::H, 0.9}, {Internal::PGV, 0.55}},
             {{Internal::Q, 0.05}},
             {{Internal::Q, 0.18076923076923068}, {Internal::H, -0.09984999999999994}}},
            {"turbine damping",
             {{External::OMEGA, 0.05}},
             {{Internal::G, 0.6},
              {Internal::Q, 0.7},
              {Internal::H, 1.1},
              {Internal::PMECH, 0.5}},
             {},
             {{Internal::PMECH, -0.2677999999999999}}},
        }};
        success *= runResidualCases(data, 0.4, turbine_cases);

        Fixture<ScalarT> curve_fixture(data);
        curve_fixture.attachAllInputs();
        success *= curve_fixture.initialize(0.33761676);
        success *= stateMatches(
            curve_fixture.hygov,
            {{Internal::C, 0.5000001394783365}, {Internal::G, 0.5000001394783365}},
            "nonidentity curve inversion");
        success *= scalarMatches(curve_fixture.input(External::PREF),
                                 0.015000004184348527,
                                 "nonidentity-curve published pref");
        success *= scalarMatches(curve_fixture.pmech(), 0.33761676, "preserved pmech value");
        success *= (curve_fixture.evaluate() == 0);
        success *= allResidualsZero(curve_fixture.hygov);

        // A flat source-curve segment must initialize to a gate on that segment.
        // makeData() uses equal power bases, At = Hdam = 1, and Qnl = 0.1,
        // so a 0.5 plateau maps to pmech = 0.4 without encoding Math::MU.
        const RealT      flat_gate_minimum = static_cast<RealT>(0.4);
        const RealT      flat_gate_maximum = static_cast<RealT>(0.6);
        const RealT      plateau_power     = static_cast<RealT>(0.5);
        const RealT      plateau_pmech     = static_cast<RealT>(0.4);
        Fixture<ScalarT> flat_fixture(makeData(),
                                      {{Params::Pgv2, plateau_power},
                                       {Params::Pgv3, plateau_power}});
        success *= flat_fixture.initialize(plateau_pmech);
        const RealT flat_gate =
            flat_fixture.hygov.y().getData()[static_cast<size_t>(Internal::G)];
        if (flat_gate < flat_gate_minimum || flat_gate > flat_gate_maximum)
        {
          std::cout << "flat-segment plateau gate "
                    << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                    << flat_gate << " is outside [" << flat_gate_minimum
                    << ", " << flat_gate_maximum << "]\n";
          success = false;
        }
        success *= (flat_fixture.evaluate() == 0);
        success *= allResidualsZero(flat_fixture.hygov);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Every Enzyme CSR row must match dependency tracking at gates inside
      /// each curve segment and at each breakpoint, and both paths must
      /// carry the PGV row's gate dependence.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto                 data = makeResidualData();
        const std::array<RealT, 9> gate_points{{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}};

        for (const RealT gate : gate_points)
        {
          const auto dependency_jacobian = dependencyTrackingJacobian(data, gate, success);
          const auto enzyme_jacobian     = enzymeJacobian(data, gate, success);

          success         *= (dependency_jacobian.size() == enzyme_jacobian.size());
          const auto rows  = std::min(dependency_jacobian.size(), enzyme_jacobian.size());
          for (size_t row = 0; row < rows; ++row)
          {
            if (!isEqual(dependency_jacobian[row], enzyme_jacobian[row], kTol))
            {
              std::cout << "HYGOV Jacobian row " << row << " at gate " << gate
                        << " mismatch between dependency tracking and Enzyme\n";
              success = false;
            }
          }

          // Guard the required PGV/G dependency even if both paths agree.
          success *= jacobianContains(
              dependency_jacobian, Internal::PGV, Internal::G, "dependency-tracking");
          success *= jacobianContains(enzyme_jacobian, Internal::PGV, Internal::G, "Enzyme");
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

      using InternalRow  = std::pair<Internal, RealT>;
      using InternalRows = std::vector<InternalRow>;
      using ExternalRow  = std::pair<External, RealT>;
      using ExternalRows = std::vector<ExternalRow>;

      /// Failure-report names for the internal rows, ordered as `Internal`.
      static constexpr std::array<const char*, static_cast<size_t>(Internal::MAXIMUM)> kRowNames{
          {"XN", "XF", "C", "G", "Q", "OMEGADB", "EF", "FC", "RC", "PGV", "H", "PMECH"}};

      struct ResidualCase
      {
        const char*  label;
        ExternalRows inputs;
        InternalRows state;
        InternalRows derivative;
        InternalRows expected;
      };

      static Data withParameters(Data                                            data,
                                 std::initializer_list<std::pair<Params, RealT>> overrides)
      {
        for (const auto& [parameter, value] : overrides)
        {
          data.parameters[parameter] = value;
        }
        return data;
      }

      /// Owns the HYGOV model, the assigned mechanical-power node, and the
      /// attached input nodes. Signal storage is declared before the model so
      /// every referenced node outlives HYGOV. Copying would invalidate the
      /// model and signal-node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, static_cast<size_t>(External::MAXIMUM)>    input_values_{};
        std::array<IdxT, static_cast<size_t>(External::MAXIMUM)> input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>,
                   static_cast<size_t>(External::MAXIMUM)>
            input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> pmech_node_;

      public:
        explicit Fixture(const Data&                                     data,
                         std::initializer_list<std::pair<Params, RealT>> overrides      = {},
                         RealT                                           system_va_base = 100.0e6)
          : hygov(withParameters(data, overrides))
        {
          hygov.setSystemBase(60.0, system_va_base);
          hygov.getSignals().template assignSignalNode<Internal::PMECH>(&pmech_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = hygov.size();

          for (size_t port = 0; port < input_values_.size(); ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = hygov.getSignals();
          signals.template attachSignalNode<External::OMEGA>(
              &input_nodes_[static_cast<size_t>(External::OMEGA)]);
          signals.template attachSignalNode<External::PREF>(
              &input_nodes_[static_cast<size_t>(External::PREF)]);
          signals.template attachSignalNode<External::PAUX>(
              &input_nodes_[static_cast<size_t>(External::PAUX)]);
        }

        /// Set the assigned mechanical-power node on the system base.
        void setPmech(RealT pmech)
        {
          pmech_node_.init(static_cast<T>(pmech));
        }

        /// Everything HYGOV initialization requires: allocation,
        /// verification, and a machine-provided mechanical-power value.
        bool prepare(RealT pmech)
        {
          const bool success = (hygov.allocate() == 0) && (hygov.verify() == 0);
          if (!success)
          {
            std::cout << "HYGOV fixture preparation failed\n";
            return false;
          }

          setPmech(pmech);
          return true;
        }

        bool initialize(RealT pmech)
        {
          if (!prepare(pmech))
          {
            return false;
          }
          if (hygov.initialize() != 0)
          {
            std::cout << "HYGOV initialization failed\n";
            return false;
          }
          return true;
        }

        int evaluate()
        {
          return hygov.evaluateResidual();
        }

        T pmech() const
        {
          return pmech_node_.read();
        }

        T& input(External port)
        {
          return input_values_[static_cast<size_t>(port)];
        }

        IdxT inputIndex(External port) const
        {
          return input_indices_[static_cast<size_t>(port)];
        }

        PhasorDynamics::Governor::Hygov<T, IdxT> hygov;
      };

      Data makeMinimalData() const
      {
        Data data;
        data.device_class          = "Hygov";
        data.disambiguation_string = "hygov_test";
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

      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(External::OMEGA) = static_cast<T>(0.02);
        fixture.input(External::PREF)  = static_cast<T>(0.31);
        fixture.input(External::PAUX)  = static_cast<T>(0.07);
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Governor::Hygov<T, IdxT>& hygov) const
      {
        setState(hygov,
                 {{Internal::XN, 0.11},
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
                  {Internal::PMECH, 0.33}});
        setDerivative(hygov,
                      {{Internal::XN, 0.01},
                       {Internal::XF, -0.02},
                       {Internal::C, 0.03},
                       {Internal::G, -0.04},
                       {Internal::Q, 0.05}});
      }

      /// Omitting every optional parameter must give exactly the model built
      /// from the defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), {}, 200.0e6);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), {}, 200.0e6);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        bool success = implicit_defaults.initialize(0.3)
                       && explicit_defaults.initialize(0.3);
        if (!success)
        {
          std::cout << "HYGOV documented-default comparison failed to initialize\n";
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
        if (!vectorUnchanged(implicit_defaults.hygov.y(),
                             copyVector(explicit_defaults.hygov.y()),
                             "documented-default state"))
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.hygov.yp(),
                             copyVector(explicit_defaults.hygov.yp()),
                             "documented-default derivative"))
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.hygov.getResidual(),
                             copyVector(explicit_defaults.hygov.getResidual()),
                             "documented-default residual"))
        {
          success = false;
        }

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.hygov);
        setAnswerKeyState(explicit_defaults.hygov);
        if (implicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (explicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (!vectorUnchanged(implicit_defaults.hygov.getResidual(),
                             copyVector(explicit_defaults.hygov.getResidual()),
                             "documented-default dynamic residual"))
        {
          success = false;
        }
        return success;
      }

      template <External variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeData());
        fixture.hygov.getSignals().template attachSignalNode<variable>(&unlinked_node);
        return fixture.hygov.verify() > 0;
      }

      template <typename VectorT>
      std::vector<RealT> copyVector(const VectorT& vector) const
      {
        const auto* values = vector.getData();
        return std::vector<RealT>(values,
                                  values + static_cast<size_t>(vector.getSize()));
      }

      template <typename VectorT>
      bool vectorUnchanged(const VectorT&            vector,
                           const std::vector<RealT>& snapshot,
                           const char*               what) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t i = 0; i < snapshot.size(); ++i)
        {
          if (!rowMatches(static_cast<RealT>(values[i]), snapshot[i], what, i, "changed"))
          {
            success = false;
          }
        }
        return success;
      }

      /// An initialization input retains exactly the value supplied by its
      /// owner, including signed infinities and NaN.
      bool scalarPreserved(RealT       actual,
                           RealT       expected,
                           const char* what,
                           size_t      row) const
      {
        bool ret = actual == expected;
        if (std::isnan(expected))
        {
          ret = std::isnan(actual);
        }
        if (!ret)
        {
          std::cout << "HYGOV " << what << " row " << row
                    << " changed mismatch: " << actual << " != " << expected << "\n";
        }
        return ret;
      }

      /// Fill the state and derivative with a recognizable ramp, then restore
      /// the aliased pmech entry, so any write by a rejected initialization
      /// is visible.
      void poisonState(Fixture<ScalarT>& fixture, RealT pmech) const
      {
        auto* y  = fixture.hygov.y().getData();
        auto* yp = fixture.hygov.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.hygov.y().getSize()); ++i)
        {
          y[i]  = 0.125 + 0.01 * static_cast<RealT>(i);
          yp[i] = -0.25 - 0.01 * static_cast<RealT>(i);
        }
        fixture.setPmech(pmech);
        fixture.hygov.y().setDataUpdated();
        fixture.hygov.yp().setDataUpdated();
      }

      /// Initialization must fail and leave the poisoned state, the seeded
      /// pmech value, and every supplied input untouched.
      bool initializationRejectedAtomically(const Data&         data,
                                            RealT               pmech,
                                            const ExternalRows& inputs,
                                            const char*         label) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        for (const auto& [port, value] : inputs)
        {
          fixture.input(port) = static_cast<ScalarT>(value);
        }
        if (!fixture.prepare(pmech))
        {
          return false;
        }

        poisonState(fixture, pmech);
        const auto y_before  = copyVector(fixture.hygov.y());
        const auto yp_before = copyVector(fixture.hygov.yp());

        bool success = true;
        if (fixture.hygov.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: " << label << "\n";
          success = false;
        }

        if (!scalarMatches(fixture.pmech(), pmech, "rejected pmech preservation"))
        {
          success = false;
        }
        for (const auto& [port, value] : inputs)
        {
          if (!scalarPreserved(static_cast<RealT>(fixture.input(port)),
                               value,
                               "external input",
                               static_cast<size_t>(port)))
          {
            success = false;
          }
        }
        if (!vectorUnchanged(fixture.hygov.y(), y_before, "state"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.hygov.yp(), yp_before, "derivative"))
        {
          success = false;
        }
        return success;
      }

      /// Write state rows and publish the update, folding in the
      /// setDataUpdated() that a hand-written write block has to remember.
      template <typename T>
      void setState(PhasorDynamics::Governor::Hygov<T, IdxT>& hygov,
                    const InternalRows&                       rows) const
      {
        auto* y = hygov.y().getData();
        for (const auto& [variable, value] : rows)
        {
          y[static_cast<size_t>(variable)] = static_cast<T>(value);
        }
        hygov.y().setDataUpdated();
      }

      template <typename T>
      void setDerivative(PhasorDynamics::Governor::Hygov<T, IdxT>& hygov,
                         const InternalRows&                       rows) const
      {
        auto* yp = hygov.yp().getData();
        for (const auto& [variable, value] : rows)
        {
          yp[static_cast<size_t>(variable)] = static_cast<T>(value);
        }
        hygov.yp().setDataUpdated();
      }

      /// Evaluate each scenario on a fresh fixture to prevent state leakage.
      template <size_t size>
      bool runResidualCases(const Data&                           data,
                            RealT                                 pmech,
                            const std::array<ResidualCase, size>& cases) const
      {
        bool success = true;
        for (const auto& test_case : cases)
        {
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          if (!fixture.initialize(pmech))
          {
            success = false;
            continue;
          }
          for (const auto& [port, value] : test_case.inputs)
          {
            fixture.input(port) = static_cast<ScalarT>(value);
          }
          setState(fixture.hygov, test_case.state);
          setDerivative(fixture.hygov, test_case.derivative);
          if (fixture.evaluate() != 0)
          {
            success = false;
          }
          if (!residualsMatch(fixture.hygov, test_case.expected, test_case.label))
          {
            success = false;
          }
        }
        return success;
      }

      /// Compare one named row and report mismatches consistently.
      bool rowMatches(RealT       actual,
                      RealT       expected,
                      const char* what,
                      size_t      row,
                      const char* context) const
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "HYGOV " << what << " row ";
        if (row < kRowNames.size())
        {
          std::cout << kRowNames[row];
        }
        else
        {
          std::cout << row;
        }
        std::cout << ' ' << context << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      template <typename VectorT, typename RowsT>
      bool rowsMatch(const VectorT& vector,
                     const RowsT&   rows,
                     const char*    what,
                     const char*    context) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (const auto& [variable, expected] : rows)
        {
          const auto row = static_cast<size_t>(variable);
          if (!rowMatches(static_cast<RealT>(values[row]), expected, what, row, context))
          {
            success = false;
          }
        }
        return success;
      }

      bool residualsMatch(const HygovT&       hygov,
                          const InternalRows& rows,
                          const char*         context = "") const
      {
        return rowsMatch(hygov.getResidual(), rows, "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const HygovT&                        hygov,
                          const std::array<InternalRow, size>& rows,
                          const char*                          context = "") const
      {
        return rowsMatch(hygov.getResidual(), rows, "residual", context);
      }

      bool stateMatches(const HygovT&       hygov,
                        const InternalRows& rows,
                        const char*         context = "") const
      {
        return rowsMatch(hygov.y(), rows, "state", context);
      }

      /// The model sits at a steady state: every residual and every
      /// derivative is zero.
      bool allResidualsZero(const HygovT& hygov) const
      {
        bool        success = true;
        const auto* f       = hygov.getResidual().getData();
        const auto* yp      = hygov.yp().getData();
        for (size_t row = 0; row < static_cast<size_t>(hygov.getResidual().getSize()); ++row)
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
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << "\n";
        return false;
      }

      void noteExpectedLogs(const char* message) const
      {
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << "\n";
        Log::setVerbosity(previous_verbosity);
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y  = fixture.hygov.y().getData();
        auto* yp = fixture.hygov.yp().getData();

        const auto model_size = static_cast<size_t>(fixture.hygov.size());
        for (size_t i = 0; i < model_size; ++i)
        {
          y[i].setVariableNumber(i);
          yp[i].setVariableNumber(i);
        }
        for (External port : {External::OMEGA, External::PREF, External::PAUX})
        {
          fixture.input(port).setVariableNumber(fixture.inputIndex(port));
        }

        fixture.hygov.y().setDataUpdated();
        fixture.hygov.yp().setDataUpdated();
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      template <typename JacobianRowsT>
      bool jacobianContains(const JacobianRowsT& rows,
                            Internal             row_variable,
                            Internal             column_variable,
                            const char*          what) const
      {
        const auto row    = static_cast<size_t>(row_variable);
        const auto column = static_cast<size_t>(column_variable);
        if (row < rows.size() && rows[row].count(column) == 1)
        {
          return true;
        }
        std::cout << "HYGOV " << what << " Jacobian row " << row
                  << " is missing column " << column << "\n";
        return false;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          const Data& data,
          RealT       gate,
          TestStatus& success) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.hygov);
        setState(fixture.hygov, {{Internal::G, gate}});
        numberVariables(fixture);
        success *= (fixture.evaluate() == 0);

        const auto                         model_size = static_cast<size_t>(fixture.hygov.size());
        std::vector<DepVar::DependencyMap> rows(model_size);
        const auto*                        f = fixture.hygov.getResidual().getData();
        for (size_t i = 0; i < model_size; ++i)
        {
          rows[i] = f[i].getDependencies();
        }
        return rows;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian(
          const Data& data,
          RealT       gate,
          TestStatus& success) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.hygov);
        setState(fixture.hygov, {{Internal::G, gate}});
        fixture.hygov.updateTime(0.0, 1.0);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.hygov.evaluateJacobian() == 0);
        success *= (fixture.hygov.constructCsr() == 0);
        return MapFromCsr(fixture.hygov.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
