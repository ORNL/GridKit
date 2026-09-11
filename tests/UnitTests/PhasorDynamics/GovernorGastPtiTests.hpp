#pragma once

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
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
        const RestoreVerbosity restore_verbosity;
        TestStatus             success = true;

        // Suppress expected errors and warnings from the invalid cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

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
        if (bound.bind(bound_y, bound_yp, bound_f, bound_abs_tol, 0) != 0)
        {
          return TestStatus(false).report(__func__);
        }
        success *= (bound.verify() == 0);

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> minimal(makeMinimalData());
        success *= (minimal.verify() > 0); // required pmech assignment is absent
        success *= (verifyData(makeData()) == 0);
        success *= defaultsMatchDocumentedValues();

        auto missing_trate = makeMinimalData();
        missing_trate.parameters.erase(Params::Trate);
        success *= (verifyData(missing_trate) > 0);

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
        Fixture<ScalarT> integer(integer_real, __func__, kTol);
        configureGastPti(integer);
        if (!integer.initialize({{Internal::PMECH, 0.4}}))
        {
          return TestStatus(false).report(__func__);
        }
        success *= integer.checkStateRows({{Internal::XFLOW, 0.8}}, "integer-valued component base");

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
          Fixture<ScalarT> invalid_base(makeData(), __func__, kTol);
          configureGastPti(invalid_base, system_base);
          if (invalid_base.model().allocate() != 0)
          {
            return TestStatus(false).report(__func__);
          }
          success *= (invalid_base.model().verify() > 0);
        }

        success *= unlinkedSignalRejected<External::OMEGA>();
        success *= unlinkedSignalRejected<External::PREF>();
        success *= aliasedSignalsRejected();

        Fixture<ScalarT> unallocated(makeData(), __func__, kTol);
        configureGastPti(unallocated);
        success *= (unallocated.model().initialize() != 0);

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

          // Keep response limits inactive while testing the time constants.
          time_data.parameters[Params::Vmin] = -40.0 / Math::MU<RealT>;
          time_data.parameters[Params::Vmax] = 1.0 + 40.0 / Math::MU<RealT>;

          Fixture<ScalarT> floors(time_data, __func__, kTol);
          configureGastPti(floors);
          if (!floors.initialize({{Internal::PMECH, 0.4}})
              || !floors.setState({
                  {Internal::XVALVE, 0.401},
                  {Internal::XFLOW, 0.4},
                  {Internal::XTEMP, 0.399},
                  {Internal::VLV, 0.402},
              })
              || !floors.evaluateResidual())
          {
            return TestStatus(false).report(__func__);
          }
          success *= floors.checkResidualRows({{Internal::XVALVE, test_case.expected_residual},
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

        Fixture<ScalarT> fixture(data, __func__, kTol);
        configureGastPti(fixture);
        attachInputs(fixture);
        // stale value the publication must replace
        if (!fixture.setInput(External::PREF, 99.0)
            || !fixture.initialize({{Internal::PMECH, 0.4}})
            || fixture.model().tagDifferentiable() != 0
            || !fixture.evaluateResidual())
        {
          return TestStatus(false).report(__func__);
        }

        const auto* y  = fixture.model().y().getData();
        success       *= scalarMatches(y[index(Internal::XVALVE)], 0.8, "XVALVE on component base");
        success       *= scalarMatches(y[index(Internal::XFLOW)], 0.8, "XFLOW on component base");
        success       *= scalarMatches(y[index(Internal::XTEMP)], 0.8, "XTEMP on component base");
        success       *= scalarMatches(y[index(Internal::VLOAD)], 0.8, "VLOAD behind the LV gate");
        success       *= scalarMatches(y[index(Internal::VTEMP)], 2.36, "VTEMP at the temperature limit");
        success       *= scalarMatches(y[index(Internal::VLV)], 0.8, "VLV at the fuel flow");
        success       *= scalarPreserved(fixture.output(Internal::PMECH), 0.4, "preserved pmech seed");

        success *= scalarPreserved(fixture.input(External::OMEGA),
                                   0.0,
                                   "preserved omega input");
        success *= scalarMatches(fixture.input(External::PREF), 0.4, "published pref");

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

        for (size_t i = 0; i < static_cast<size_t>(fixture.model().size()); ++i)
        {
          const bool expected = i <= index(Internal::XTEMP);
          if (fixture.model().tag()[i] != expected)
          {
            std::cout << "GASTPTI differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= fixture.checkSteadyState();

        // A system-base reference step lands on the droop row scaled by the
        // base ratio.
        // the published 0.4 plus a 0.1 step
        if (!fixture.setInput(External::PREF, 0.5)
            || !fixture.evaluateResidual())
        {
          return TestStatus(false).report(__func__);
        }
        success *= fixture.checkResidualRows({{Internal::VLOAD, 0.01}}, "reference step on the component base");

        // GridKit deliberately leaves references above At uncapped.
        // 2.2 on component base; At = 2.0
        if (!fixture.setInput(External::PREF, 1.1)
            || !fixture.evaluateResidual())
        {
          return TestStatus(false).report(__func__);
        }
        success *= fixture.checkResidualRows({{Internal::VLOAD, 0.07}}, "uncapped reference above At");

        // Unattached ports fall back to the reference latched by
        // initialize(), so the same steady state holds without a controller.
        Fixture<ScalarT> latched(data, __func__, kTol);
        configureGastPti(latched);
        if (!latched.initialize({{Internal::PMECH, 0.4}})
            || !latched.evaluateResidual())
        {
          return TestStatus(false).report(__func__);
        }
        success *= latched.checkSteadyState();

        constexpr RealT initial_pmech       = 0.4;
        constexpr RealT system_to_component = 2.0;
        constexpr RealT component_to_system = ONE<RealT> / system_to_component;
        constexpr RealT droop               = 0.05;
        constexpr RealT temperature_limit   = 2.0;
        constexpr RealT temperature_gain    = 0.3;
        constexpr RealT turbine_damping     = 0.1;

        const std::array<RealT, 2> speed_cases{{0.05, -0.05}};
        for (const RealT omega : speed_cases)
        {
          const RealT xflow = system_to_component * initial_pmech
                              + turbine_damping * omega;
          const RealT vtemp = temperature_limit
                              + temperature_gain * (temperature_limit - xflow);

          Fixture<ScalarT> speed_fixture(data, __func__, kTol);
          configureGastPti(speed_fixture);
          attachInputs(speed_fixture);
          if (!speed_fixture.setInput(External::OMEGA, omega)
              || !speed_fixture.initialize({{Internal::PMECH, initial_pmech}}))
          {
            return TestStatus(false).report(__func__);
          }
          success             *= speed_fixture.checkStateRows({
                                                      {Internal::XVALVE, xflow},
                                                      {Internal::XFLOW, xflow},
                                                      {Internal::XTEMP, xflow},
                                                      {Internal::VTEMP, vtemp},
                                                      {Internal::VLV, xflow},
                                                  },
                                                  "signed nonzero-speed initialization");
          success             *= scalarPreserved(speed_fixture.output(Internal::PMECH),
                                     initial_pmech,
                                     "signed-speed pmech preservation");
          success             *= scalarPreserved(speed_fixture.input(External::OMEGA),
                                     omega,
                                     "signed-speed input preservation");
          const auto* speed_y  = speed_fixture.model().y().getData();
          const RealT vload    = static_cast<RealT>(speed_y[index(Internal::VLOAD)]);
          const RealT pref     = component_to_system * (vload + omega / droop);
          success             *= scalarMatches(speed_fixture.input(External::PREF),
                                   pref,
                                   "signed-speed pref publication");
          if (!speed_fixture.evaluateResidual())
          {
            return TestStatus(false).report(__func__);
          }
          success *= speed_fixture.checkSteadyState();
        }

        return success.report(__func__);
      }

      /// Initialization-domain boundaries, effective limits, and exact failure
      /// atomicity.
      TestOutcome initializationDomain()
      {
        const RestoreVerbosity restore_verbosity;
        TestStatus             success = true;

        // Suppress expected errors and response-limit warnings from the cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

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
        Fixture<ScalarT> invalid(invalid_data, __func__, kTol);
        configureGastPti(invalid);
        attachInputs(invalid);
        if (invalid.prepare() || !poisonState(invalid, 0.4))
        {
          return TestStatus(false).report(__func__);
        }
        const auto invalid_before = invalid.snapshot();
        if (invalid.model().initialize() == 0)
        {
          std::cout << "Expected initialization rejection: invalid configuration\n";
          success = false;
        }
        success *= invalid.checkUnchanged(invalid_before);

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
        Fixture<ScalarT> over_rated(makeResidualData(), __func__, kTol);
        // fuel flow 1.2 above Vmax = 1.1
        if (!initializeGastPti(over_rated, 0.6))
        {
          return TestStatus(false).report(__func__);
        }
        success *= over_rated.checkStateRows({{Internal::XVALVE, 1.2},
                                              {Internal::XFLOW, 1.2},
                                              {Internal::VLV, 1.2}},
                                             "over-rated dispatch");
        success *= scalarPreserved(over_rated.output(Internal::PMECH),
                                   0.6,
                                   "preserved over-rated pmech seed");
        success *= over_rated.checkSteadyState();

        // A failed reinitialization must preserve the last committed effective
        // limits as well as state, derivatives, and pref.
        constexpr RealT over_rated_pmech    = 0.6;
        constexpr RealT valve_time_constant = 0.35;
        constexpr RealT boundary_weight     = 0.5;
        const RealT     boundary_command    = 40.0 / Math::MU<RealT>;

        auto data_reused                   = makeResidualData();
        data_reused.parameters[Params::T1] = valve_time_constant;
        Fixture<ScalarT> reused(data_reused, __func__, kTol);
        if (!initializeGastPti(reused, over_rated_pmech)
            || !reused.setState({{Internal::PMECH, 1.0}}))
        {
          return TestStatus(false).report(__func__);
        }
        const auto  before         = reused.snapshot();
        const RealT upper_boundary = reused.state(Internal::XVALVE);
        if (reused.model().initialize() == 0)
        {
          std::cout << "Expected failed GASTPTI reinitialization\n";
          success = false;
        }
        success *= reused.checkUnchanged(before);
        if (!reused.setState({{Internal::PMECH, over_rated_pmech}}))
        {
          return TestStatus(false).report(__func__);
        }
        if (!reused.setState({{Internal::XVALVE, upper_boundary},
                              {Internal::VLV, upper_boundary + boundary_command}})
            || !reused.setDerivative({{Internal::XVALVE, 0.0}})
            || !reused.evaluateResidual())
        {
          return TestStatus(false).report(__func__);
        }
        const RealT expected_boundary_response =
            boundary_weight * boundary_command / valve_time_constant;
        success *= reused.checkResidualRows({{Internal::XVALVE, expected_boundary_response}}, "failed reinitialization preserves effective limits");

        struct EqualLimitTemperatureCase
        {
          const char* label;
          RealT       at;
          RealT       vtemp;
        };

        const std::array<EqualLimitTemperatureCase, 2> equal_limit_temperature_cases{{
            {"equal valve limits at the temperature limit", 0.8, 0.8},
            {"equal valve limits above the temperature limit", 0.0, -0.32},
        }};
        for (const auto& test_case : equal_limit_temperature_cases)
        {
          auto equal_limit_data                     = makeResidualData();
          equal_limit_data.parameters[Params::At]   = test_case.at;
          equal_limit_data.parameters[Params::Vmin] = 0.8;
          equal_limit_data.parameters[Params::Vmax] = 0.8;

          Fixture<ScalarT> equal_limits(equal_limit_data, __func__, kTol);
          if (!initializeGastPti(equal_limits, 0.4))
          {
            return TestStatus(false).report(__func__);
          }
          success *= equal_limits.checkStateRows({
                                                     {Internal::XVALVE, 0.8},
                                                     {Internal::XFLOW, 0.8},
                                                     {Internal::XTEMP, 0.8},
                                                     {Internal::VLOAD, 0.8},
                                                     {Internal::VTEMP, test_case.vtemp},
                                                     {Internal::PMECH, 0.4},
                                                 },
                                                 test_case.label);
          if (!equal_limits.evaluateResidual())
          {
            return TestStatus(false).report(__func__);
          }
          success *= equal_limits.checkSteadyState();
        }

        // An unattached reference retains its last successful latch when a
        // later active reinitialization is rejected.
        Fixture<ScalarT> latched(makeResidualData(), __func__, kTol);
        configureGastPti(latched);
        if (!latched.initialize({{Internal::PMECH, 0.4}})
            || !latched.setState({{Internal::PMECH, 1.0}}))
        {
          return TestStatus(false).report(__func__);
        }
        success *= (latched.model().initialize() != 0);
        if (!latched.setState({{Internal::PMECH, 0.4}})
            || !latched.evaluateResidual())
        {
          return TestStatus(false).report(__func__);
        }
        success *= latched.checkSteadyState();

        // A zero mechanical-power seed stays admissible.
        Fixture<ScalarT> zero_seed(makeResidualData(), __func__, kTol);
        if (!initializeGastPti(zero_seed, 0.0))
        {
          return TestStatus(false).report(__func__);
        }
        success *= zero_seed.checkStateRows({{Internal::XFLOW, 0.0}, {Internal::VTEMP, 2.52}}, "zero seed");
        success *= zero_seed.checkSteadyState();

        auto negative_data                     = makeResidualData();
        negative_data.parameters[Params::Vmin] = -1.0;
        Fixture<ScalarT> negative_seed(negative_data, __func__, kTol);
        if (!initializeGastPti(negative_seed, -0.1))
        {
          return TestStatus(false).report(__func__);
        }
        success *= negative_seed.checkStateRows({{Internal::XFLOW, -0.2},
                                                 {Internal::VLOAD, -0.2},
                                                 {Internal::PMECH, -0.1}},
                                                "negative finite dispatch");
        success *= negative_seed.checkSteadyState();

        return success.report(__func__);
      }

      /// A fixed near-closed temperature-gate case proves that initialization
      /// uses the inverse smooth ramp and rests all seven residuals exactly
      /// within the documented behavior tolerance.
      TestOutcome initializationExactness()
      {
        TestStatus success = true;

        constexpr RealT initial_flow        = 0.8;
        constexpr RealT system_to_component = 2.0;
        constexpr RealT initial_pmech       = initial_flow / system_to_component;
        constexpr RealT temperature_gain    = 0.4;
        const RealT     temperature_margin  = 0.02 / Math::MU<RealT>;

        auto data                   = makeResidualData();
        data.parameters[Params::Kt] = temperature_gain;
        data.parameters[Params::At] =
            initial_flow + temperature_margin / (ONE<RealT> + temperature_gain);

        Fixture<ScalarT> fixture(data, __func__, kTol);
        if (!initializeGastPti(fixture, initial_pmech))
        {
          return TestStatus(false).report(__func__);
        }
        success *= fixture.checkStateRows({{Internal::VTEMP, initial_flow + temperature_margin},
                                           {Internal::VLV, initial_flow}},
                                          "near-gate initialization");

        const auto* y     = fixture.model().y().getData();
        const RealT vload = static_cast<RealT>(y[index(Internal::VLOAD)]);
        const RealT vtemp = static_cast<RealT>(y[index(Internal::VTEMP)]);
        if (!(vload > vtemp))
        {
          std::cout << "GASTPTI near-gate initialization selected the wrong demand side\n";
          success = false;
        }
        success *= fixture.checkSteadyState();

        // A very large but finite temperature margin must not erase the
        // ordinary-sized load demand through catastrophic cancellation.
        // With Kt = 0, the exact initialized load demand remains xF0 = 0.8.
        auto large_margin_data                   = makeResidualData();
        large_margin_data.parameters[Params::At] = 1.0e16;
        large_margin_data.parameters[Params::Kt] = 0.0;
        Fixture<ScalarT> large_margin(large_margin_data, __func__, kTol);
        if (!initializeGastPti(large_margin, 0.4))
        {
          return TestStatus(false).report(__func__);
        }
        success *= large_margin.checkStateRows({{Internal::VLOAD, 0.8},
                                                {Internal::VLV, 0.8}},
                                               "large finite temperature margin");
        success *= large_margin.checkSteadyState();

        return success.report(__func__);
      }

      /// Check all seven equations against arithmetic and ideal limiter values.
      TestOutcome residualEquations()
      {
        Fixture<ScalarT> fixture(makeResidualData(), __func__, kTol);
        if (!initializeGastPti(fixture, 0.4) || !fixture.setPoint(residualPoint()))
          return TestStatus(false).report(__func__);
        TestStatus success = fixture.checkResidualRows({
            {Internal::XFLOW, 0.22},
            {Internal::XTEMP, 0.07},
            {Internal::VLOAD, -0.0326},
            {Internal::VTEMP, 1.0},
            {Internal::PMECH, -0.1424},
        });

        success *= scalarMatches(fixture.residual(Internal::XVALVE), 0.19, "interior valve rate", kTol + 0.6 * std::exp(-0.49 * Math::MU<RealT>));
        success *= scalarMatches(fixture.residual(Internal::VLV), 0.15, "load demand", kTol + std::exp(-0.57 * Math::MU<RealT>) / Math::MU<RealT>);
        return success.report(__func__);
      }

      /// Valve anti-windup, speed/damping signs, and adjusted Normal limits.
      TestOutcome governorControl()
      {
        const RestoreVerbosity restore_verbosity;
        TestStatus             success = true;

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
            {"Vmax admits a restoring valve rate", 1.6, 1.25, -1.0},
            {"Vmin admits a restoring valve rate", -0.45, -0.1, 1.0},
        }};
        for (const auto& test_case : antiwindup_cases)
        {
          Fixture<ScalarT> antiwindup(makeResidualData(), __func__, kTol + 3.0 * std::exp(-0.25 * Math::MU<RealT>));

          success *= initializeGastPti(antiwindup, 0.4)
                     && antiwindup.checkResidualRows(
                         {.state = {{Internal::XVALVE, test_case.xvalve}, {Internal::VLV, test_case.vlv}}, .derivative = {{Internal::XVALVE, 0.0}}},
                         {{Internal::XVALVE, test_case.expected}},
                         test_case.label);
        }

        // A speed deviation enters the droop and turbine-damping rows.
        Fixture<ScalarT> speed_step(makeResidualData(), __func__, kTol);

        success *= initializeGastPti(speed_step, 0.4)
                   && speed_step.checkResidualRows(
                       {.inputs = {{External::OMEGA, 0.05}}},
                       {{Internal::VLOAD, -0.05},
                        {Internal::PMECH, -0.006}},
                       "speed deviation in the droop and damping rows");

        // Normal response expands both sides of the configured interval to
        // admit the initialized flow. The derived boundary must be used thereafter.
        constexpr RealT valve_time_constant = 0.35;
        constexpr RealT boundary_weight     = 0.5;
        const RealT     command_magnitude   = 40.0 / Math::MU<RealT>;
        constexpr RealT over_rated_pmech    = 0.6;

        auto response_data                   = makeResidualData();
        response_data.parameters[Params::T1] = valve_time_constant;

        struct EffectiveBoundaryCase
        {
          const char* label;
          RealT       pmech;
          RealT       command;
        };

        const std::array<EffectiveBoundaryCase, 2> effective_boundary_cases{{
            {"adjusted upper response boundary", over_rated_pmech, command_magnitude},
            {"adjusted lower response boundary", ZERO<RealT>, -command_magnitude},
        }};
        // Suppress expected response-limit adjustment warnings from these cases.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);
        for (const auto& test_case : effective_boundary_cases)
        {
          Fixture<ScalarT> response(response_data, __func__, kTol);
          if (!initializeGastPti(response, test_case.pmech))
          {
            return TestStatus(false).report(__func__);
          }

          const RealT boundary = response.state(Internal::XVALVE);
          const RealT expected =
              boundary_weight * test_case.command / valve_time_constant;
          success *= response.checkResidualRows(
              {.state      = {{Internal::XVALVE, boundary}, {Internal::VLV, boundary + test_case.command}},
               .derivative = {{Internal::XVALVE, ZERO<RealT>}}},
              {{Internal::XVALVE, expected}},
              test_case.label);
        }

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
            {"equal demands split the smooth LV gate", 0.9, 0.9, 0.9 - std::log(2.0) / Math::MU<RealT>},
        }};
        for (const auto& test_case : gate_cases)
        {
          // Away from equality, the soft-min error is bounded by its exponential tail.
          const RealT      gap       = std::abs(test_case.vload - test_case.vtemp);
          const RealT      tolerance = gap == 0.0 ? kTol
                                                  : kTol + std::exp(-Math::MU<RealT> * gap) / Math::MU<RealT>;
          Fixture<ScalarT> gate(makeResidualData(), __func__, tolerance);

          success *= initializeGastPti(gate, 0.4)
                     && gate.checkResidualRows(
                         {.state = {{Internal::VLOAD, test_case.vload},
                                    {Internal::VTEMP, test_case.vtemp},
                                    {Internal::VLV, 0.0}}},
                         {{Internal::VLV, test_case.expected}},
                         test_case.label);
        }

        // The exhaust-temperature feedback drives the temperature demand.
        Fixture<ScalarT> feedback(makeResidualData(), __func__, kTol);

        success *= initializeGastPti(feedback, 0.4)
                   && feedback.checkResidualRows(
                       {.state = {{Internal::XTEMP, 0.9}, {Internal::VTEMP, 1.1}}},
                       {{Internal::VTEMP, 1.06}},
                       "temperature feedback");

        // At equality, the smooth low-value selector splits its sensitivity
        // evenly between the two demand signals.
        Fixture<DependencyTracking::Variable> selector(makeResidualData(), __func__, kTol);
        if (!initializeGastPti(selector, 0.4)
            || !selector.setState({{Internal::VLOAD, 0.9},
                                   {Internal::VTEMP, 0.9},
                                   {Internal::VLV, 0.7}}))
        {
          return TestStatus(false).report(__func__);
        }

        success *= selector.checkJacobianRow(
            Internal::VLV,
            {{Internal::VLOAD, 0.5}, {Internal::VTEMP, 0.5}, {Internal::VLV, -1.0}});

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Enzyme and dependency tracking agree row-for-row across selector,
      /// anti-windup, and collapsed-limit configurations.
      TestOutcome jacobian()
      {
        const RestoreVerbosity restore_verbosity;
        TestStatus             success = true;

        constexpr RealT initial_pmech    = 0.4;
        constexpr RealT over_rated_pmech = 0.6;
        constexpr RealT boundary_command = 0.25;
        constexpr RealT collapsed_limit  = 0.8;

        const auto data = makeResidualData();

        const auto compare = [&](const Data&   case_data,
                                 RealT         pmech,
                                 const char*   context,
                                 const Values& overrides,
                                 bool          retry = false)
        {
          return Fixture<ScalarT>::checkJacobian(
              case_data,
              [&](auto& fixture)
              {
                const RealT seed = retry ? initial_pmech : pmech;
                if (!initializeGastPti(fixture, seed))
                  return false;
                if (retry
                    && (!fixture.setInput(External::PREF, 0.0)
                        || !fixture.initialize({{Internal::PMECH, pmech}})))
                  return false;
                return fixture.setPoint(residualPoint()) && fixture.setState(overrides);
              },
              {0.0, 1.0, 2.5},
              context,
              kTol);
        };

        success *= compare(data,
                           initial_pmech,
                           "load-limited Enzyme versus dependency tracking",
                           {});
        success *= compare(data,
                           initial_pmech,
                           "temperature-limited Enzyme versus dependency tracking",
                           {{Internal::VLOAD, 1.5}, {Internal::VTEMP, 0.3}});
        success *= compare(data,
                           initial_pmech,
                           "equal-selector Enzyme versus dependency tracking",
                           {{Internal::VLOAD, 0.9}, {Internal::VTEMP, 0.9}});
        success *= compare(data,
                           initial_pmech,
                           "blocked-response Enzyme versus dependency tracking",
                           {{Internal::XVALVE, 1.6}, {Internal::VLV, 1.85}});
        success *= compare(data,
                           initial_pmech,
                           "restoring-response Enzyme versus dependency tracking",
                           {{Internal::XVALVE, 1.6}, {Internal::VLV, 1.35}});

        // Suppress the expected response-limit adjustment warning for this case.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        Fixture<ScalarT> adjusted(data, __func__, kTol);
        if (!initializeGastPti(adjusted, over_rated_pmech))
        {
          return TestStatus(false).report(__func__);
        }

        const RealT adjusted_boundary  = adjusted.state(Internal::XVALVE);
        success                       *= compare(data,
                           over_rated_pmech,
                           "adjusted-boundary Enzyme versus dependency tracking",
                                                 {{Internal::XVALVE, adjusted_boundary},
                                                  {Internal::VLV, adjusted_boundary + boundary_command}});

        auto collapsed_data                     = data;
        collapsed_data.parameters[Params::Vmin] = collapsed_limit;
        collapsed_data.parameters[Params::Vmax] = collapsed_limit;

        // Suppress the expected response-limit adjustment warning for this case.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        success *= compare(collapsed_data,
                           initial_pmech,
                           "collapsed Enzyme versus dependency tracking",
                           {});
        success *= compare(collapsed_data,
                           over_rated_pmech,
                           "reinitialized Enzyme versus dependency tracking",
                           {},
                           true);

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

      template <typename T>
      using Fixture = ComponentTestFixture<PhasorDynamics::Governor::GastPti, T, IdxT>;

      using Values = typename Fixture<ScalarT>::Values;
      using Point  = typename Fixture<ScalarT>::Point;

      /// Model-specific bases and required output; optional inputs stay explicit.
      template <typename T>
      void configureGastPti(Fixture<T>& fixture, RealT system_va_base = 100.0e6) const
      {
        fixture.model().setSystemBase(60.0, system_va_base);
        fixture.template assignOutput<Internal::PMECH>();
      }

      template <typename T>
      void attachInputs(Fixture<T>& fixture) const
      {
        fixture.template attachInput<External::OMEGA>(0.0);
        fixture.template attachInput<External::PREF>(0.0);
      }

      template <typename T>
      bool initializeGastPti(Fixture<T>& fixture, RealT pmech) const
      {
        configureGastPti(fixture);
        attachInputs(fixture);
        return fixture.initialize({{Internal::PMECH, pmech}});
      }

      /// Restore expected-error log suppression on prerequisite failure too.
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
        data.device_class              = "GastPti";
        data.disambiguation_string     = "gastpti_test";
        data.parameters[Params::Trate] = 100.0;
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

      /// A common point for both scalar types, with distinct values in every row.
      Point residualPoint() const
      {
        return {
            .inputs = {{External::OMEGA, 0.02}, {External::PREF, 0.31}},
            .state  = {
                {Internal::XVALVE, 0.61},
                {Internal::XFLOW, 0.52},
                {Internal::XTEMP, 0.3},
                {Internal::VLOAD, 0.83},
                {Internal::VTEMP, 1.4},
                {Internal::VLV, 0.68},
                {Internal::PMECH, 0.33},
            },
            .derivative = {{Internal::XVALVE, 0.01}, {Internal::XFLOW, -0.02}, {Internal::XTEMP, 0.03}}};
      }

      /// Omitting every optional parameter must give exactly the model built
      /// from the defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), __func__, kTol);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), __func__, kTol);
        if (!initializeGastPti(implicit_defaults, 0.3)
            || !initializeGastPti(explicit_defaults, 0.3)
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
          return false;
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
          return false;
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
          return false;
        }
        if (!(input_alias.verify() > 0))
        {
          success = false;
        }

        return success;
      }

      /// Fill the state and derivative with a recognizable ramp, then re-seed
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

      bool initializationRejectedAtomically(const Data& data,
                                            RealT       pmech,
                                            const char* label,
                                            RealT       omega = ZERO<RealT>) const
      {
        Fixture<ScalarT> fixture(data, __func__, kTol);
        configureGastPti(fixture);
        attachInputs(fixture);
        // The reference must stay untouched on rejection.
        if (!fixture.setInput(External::OMEGA, omega)
            || !fixture.setInput(External::PREF, 77.0)
            || !fixture.prepare()
            || !poisonState(fixture, pmech))
        {
          return false;
        }

        const auto before   = fixture.snapshot();
        const bool rejected = fixture.model().initialize() != 0;
        if (!rejected)
          std::cout << "Expected initialization rejection: " << label << '\n';
        const bool unchanged = fixture.checkUnchanged(before);
        const bool success   = rejected && unchanged;
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
    };
  } // namespace Testing
} // namespace GridKit
