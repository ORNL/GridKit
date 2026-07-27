#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
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

      // GASTPTI initialization seats the load demand behind the smooth LV
      // gate through the exact inverse smooth ramp, so steady residuals rest
      // at O(1e-15). Behavioral comparisons use a tolerance well above that
      // gap and well below every pinned answer-key digit.
      static constexpr RealT kBehaviorTol = 1.0e-9;

      // Enzyme and dependency tracking traverse the same smooth expressions
      // differently; their double-precision derivatives agree to O(1e-10).
      static constexpr RealT kJacobianTol = 1.0e-9;

      /// Construction and every verify() error class, including parameter
      /// types, parameter relationships, the mode-dependent valve limits,
      /// and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> empty;
        success *= (empty.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() == 0);

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> configured(makeData());
        success *= (configured.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() == 0);

        noteExpectedLogs("Testing GASTPTI defaults and invalid configurations. "
                         "Logged errors and mode/limit/time-constant warnings are expected.");

        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> minimal(makeMinimalData());
        success *= (minimal.verify() == 0);
        success *= defaultsMatchDocumentedValues();

        success *= invalidParameterCase(Params::R, 0.0);
        success *= invalidParameterCase(Params::T1, -0.1);
        success *= invalidParameterCase(Params::T2, -0.1);
        success *= invalidParameterCase(Params::T3, -0.1);
        success *= invalidParameterCase(Params::At, -0.1);
        success *= invalidParameterCase(Params::Kt, -0.1);
        success *= invalidParameterCase(Params::Vmin, 1.2); // equals Vmax in Normal mode
        success *= invalidParameterCase(Params::Vmin, 2.0); // above Vmax
        success *= invalidParameterCase(Params::Dturb, -0.1);
        success *= invalidParameterCase(Params::Trate, 0.0);

        // The response mode must be an integer inside {0, 1, 2}.
        auto real_mode                     = makeData();
        real_mode.parameters[Params::mode] = 2.0;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> real_mode_model(real_mode);
        success *= (real_mode_model.verify() > 0);

        auto invalid_mode                     = makeData();
        invalid_mode.parameters[Params::mode] = static_cast<IdxT>(3);
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> invalid_mode_model(invalid_mode);
        success *= (invalid_mode_model.verify() > 0);

        // Down Only is accepted with a warning and simulated as Normal.
        auto down_only                     = makeData();
        down_only.parameters[Params::mode] = static_cast<IdxT>(Mode::DownOnly);
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> down_only_model(down_only);
        success *= (down_only_model.verify() == 0);

        // Fixed mode admits equal valve limits; Normal mode rejected them above.
        auto fixed_equal                     = makeData();
        fixed_equal.parameters[Params::mode] = static_cast<IdxT>(Mode::Fixed);
        fixed_equal.parameters[Params::Vmin] = 0.5;
        fixed_equal.parameters[Params::Vmax] = 0.5;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> fixed_equal_model(fixed_equal);
        success *= (fixed_equal_model.verify() == 0);

        // Narrow valve limits warn without failing verification.
        auto narrow                     = makeData();
        narrow.parameters[Params::Vmax] = 0.01;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> narrow_model(narrow);
        success *= (narrow_model.verify() == 0);

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                      = makeData();
        integer_real.parameters[Params::Trate] = static_cast<IdxT>(100);
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> integer_real_model(integer_real);
        success *= (integer_real_model.verify() == 0);

        auto bad_numeric_type                   = makeData();
        bad_numeric_type.parameters[Params::T1] = true;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> bad_numeric_model(bad_numeric_type);
        success *= (bad_numeric_model.verify() > 0);

        success *= unlinkedSignalRejected<Ext::OMEGA>();
        success *= unlinkedSignalRejected<Ext::PREF>();

        // All three zero time constants use the documented numerical floor
        // and still admit a consistent steady-state initialization.
        auto zero_time                   = makeData();
        zero_time.parameters[Params::T1] = 0.0;
        zero_time.parameters[Params::T2] = 0.0;
        zero_time.parameters[Params::T3] = 0.0;

        Fixture<ScalarT> fixture(zero_time);
        success *= fixture.initialize(0.4);
        success *= (fixture.evaluate() == 0);
        success *= allResidualsZero(fixture.gastpti);

        return success.report(__func__);
      }

      /// A nonidentity power-base initialization with every port attached.
      /// The machine-seeded pmech node must remain unchanged while GASTPTI
      /// initializes and publishes its resolved load reference.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        auto data                      = makeData();
        data.parameters[Params::Trate] = 50.0;

        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        fixture.input(E::PREF)  = 99.0; // stale value the publication must replace
        success                *= fixture.initialize(0.4);
        success                *= (fixture.gastpti.tagDifferentiable() == 0);
        success                *= (fixture.evaluate() == 0);

        const auto* y  = fixture.gastpti.y().getData();
        success       *= scalarMatches(y[I::XVALVE], 0.8, "XVALVE on component base");
        success       *= scalarMatches(y[I::XFLOW], 0.8, "XFLOW on component base");
        success       *= scalarMatches(y[I::XTEMP], 0.8, "XTEMP on component base");
        success       *= scalarMatches(y[I::VLOAD], 0.8, "VLOAD behind the LV gate");
        success       *= scalarMatches(y[I::VTEMP], 2.36, "VTEMP at the temperature limit");
        success       *= scalarMatches(y[I::VLV], 0.8, "VLV at the fuel flow");
        success       *= scalarMatches(fixture.pmech(), 0.4, "preserved pmech seed");

        success *= scalarMatches(fixture.input(E::OMEGA), 0.0, "preserved omega input");
        success *= scalarMatches(fixture.input(E::PREF), 0.4, "published pref");

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
          const bool expected = i <= I::XTEMP;
          if (fixture.gastpti.tag()[i] != expected)
          {
            std::cout << "GASTPTI differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= allResidualsZero(fixture.gastpti);

        // A system-base reference step lands on the droop row scaled by the
        // base ratio.
        fixture.input(E::PREF)  = 0.5; // the published 0.4 plus a 0.1 step
        success                *= (fixture.evaluate() == 0);
        success                *= residualsMatch(fixture.gastpti,
                                                 {{I::VLOAD, 0.01}},
                                  "reference step on the component base");

        // Unattached ports fall back to the reference latched by
        // initialize(), so the same steady state holds without a controller.
        Fixture<ScalarT> fallback(data);
        success *= fallback.initialize(0.4);
        success *= (fallback.evaluate() == 0);
        success *= allResidualsZero(fallback.gastpti);

        return success.report(__func__);
      }

      /// Temperature-gate and configuration initialization domains. Every
      /// rejection is atomic; an over-rated dispatch, a Fixed-mode point
      /// outside its limits, and a zero power seed remain admissible.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible GASTPTI temperature-gate and "
                         "configuration initialization points. Logged errors "
                         "and rating warnings are expected.");

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

        // An over-rated Normal-mode dispatch warns but initializes at the
        // dispatched value, resting against the saturated valve limit.
        Fixture<ScalarT> over_rated(makeResidualData());
        over_rated.attachAllInputs();
        success *= over_rated.initialize(0.6); // fuel flow 1.2 above Vmax = 1.1
        success *= stateMatches(over_rated.gastpti,
                                {{I::XVALVE, 1.2}, {I::XFLOW, 1.2}, {I::VLV, 1.2}},
                                "over-rated dispatch");
        success *= scalarMatches(over_rated.pmech(), 0.6, "preserved over-rated pmech seed");
        success *= (over_rated.evaluate() == 0);
        success *= allResidualsZero(over_rated.gastpti);

        // A Fixed-mode point outside its equal valve limits holds the
        // dispatch with frozen turbine dynamics.
        auto fixed_data                     = makeResidualData();
        fixed_data.parameters[Params::mode] = static_cast<IdxT>(Mode::Fixed);
        fixed_data.parameters[Params::Vmin] = 0.3;
        fixed_data.parameters[Params::Vmax] = 0.3;
        Fixture<ScalarT> fixed_fixture(fixed_data);
        fixed_fixture.attachAllInputs();
        success *= fixed_fixture.initialize(0.4);
        success *= (fixed_fixture.evaluate() == 0);
        success *= allResidualsZero(fixed_fixture.gastpti);

        // A zero mechanical-power seed stays admissible.
        Fixture<ScalarT> zero_seed(makeResidualData());
        zero_seed.attachAllInputs();
        success *= zero_seed.initialize(0.0);
        success *= stateMatches(zero_seed.gastpti,
                                {{I::XFLOW, 0.0}, {I::VTEMP, 2.52}},
                                "zero seed");
        success *= (zero_seed.evaluate() == 0);
        success *= allResidualsZero(zero_seed.gastpti);

        return success.report(__func__);
      }

      /// A fixed numerical answer key for all 7 GASTPTI residual rows. The
      /// expected values are literals, not a second implementation of
      /// GASTPTI.
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
        const std::array<Row, I::MAXIMUM> expected{{
            {I::XVALVE, 0.24614285714285705},
            {I::XFLOW, 0.17755555555555544},
            {I::XTEMP, -0.001181818181818159},
            {I::VLOAD, -0.0326},
            {I::VTEMP, 0.978},
            {I::VLV, 0.12},
            {I::PMECH, -0.11239999999999999},
        }};

        success *= (static_cast<size_t>(fixture.gastpti.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.gastpti, expected);

        return success.report(__func__);
      }

      /// Valve anti-windup at every controller direction, speed deviation in
      /// the droop and damping rows, and the Fixed and Down Only response
      /// modes against the Normal-mode answer key.
      TestOutcome governorControl()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);

        // The valve anti-windup at three controller directions: both
        // saturations block an outward rate and Vmax admits a restoring one.
        struct AntiWindupCase
        {
          const char* label;
          RealT       xvalve;
          RealT       vlv;
          RealT       expected;
        };

        for (const auto& test_case : std::array<AntiWindupCase, 3>{{
                 {"Vmax blocks an outward valve rate", 1.6, 1.85, 0.0},
                 {"Vmin blocks an outward valve rate", -0.45, -0.7, 0.0},
                 {"Vmax admits a restoring valve rate", 1.6, 1.35, -0.7142857142857143},
             }})
        {
          setState(fixture.gastpti,
                   {{I::XVALVE, test_case.xvalve}, {I::VLV, test_case.vlv}});
          setDerivative(fixture.gastpti, {{I::XVALVE, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.gastpti,
                                    {{I::XVALVE, test_case.expected}},
                                    test_case.label);
        }

        // A speed deviation enters the droop and turbine-damping rows.
        fixture.input(E::OMEGA)  = 0.05;
        success                 *= (fixture.evaluate() == 0);
        success                 *= residualsMatch(fixture.gastpti,
                                                  {{I::VLOAD, -0.05}, {I::PMECH, -0.006}},
                                  "speed deviation in the droop and damping rows");
        fixture.input(E::OMEGA)  = 0.0;

        // Fixed mode freezes the three turbine states at the answer-key
        // point while the algebraic rows keep their Normal-mode values.
        auto fixed_data                     = makeResidualData();
        fixed_data.parameters[Params::mode] = static_cast<IdxT>(Mode::Fixed);
        Fixture<ScalarT> fixed_fixture(fixed_data);
        fixed_fixture.attachAllInputs();
        success *= fixed_fixture.initialize(0.4);
        setAnswerKeyInputs(fixed_fixture);
        setAnswerKeyState(fixed_fixture.gastpti);
        success *= (fixed_fixture.evaluate() == 0);
        success *= residualsMatch(fixed_fixture.gastpti,
                                  {{I::XVALVE, -0.011}, {I::XFLOW, 0.022}, {I::XTEMP, -0.033}},
                                  "Fixed mode freezes the turbine states");

        // Down Only warns and reproduces the Normal-mode answer key.
        noteExpectedLogs("Testing the GASTPTI Down Only response mode. "
                         "The logged mode warning is expected.");
        auto down_only_data                     = makeResidualData();
        down_only_data.parameters[Params::mode] = static_cast<IdxT>(Mode::DownOnly);
        Fixture<ScalarT> down_only_fixture(down_only_data);
        down_only_fixture.attachAllInputs();
        success *= down_only_fixture.initialize(0.4);
        setAnswerKeyInputs(down_only_fixture);
        setAnswerKeyState(down_only_fixture.gastpti);
        success *= (down_only_fixture.evaluate() == 0);
        success *= residualsMatch(down_only_fixture.gastpti,
                                  {{I::XVALVE, 0.24614285714285705},
                                   {I::XFLOW, 0.17755555555555544},
                                   {I::XTEMP, -0.001181818181818159}},
                                  "Down Only simulated as Normal");

        return success.report(__func__);
      }

      /// The smooth LV gate on both demand sides and at demand equality, the
      /// exhaust-temperature feedback row, and initialization against a
      /// near-closed temperature gate.
      TestOutcome temperatureLimiting()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);

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
          setState(fixture.gastpti,
                   {{I::VLOAD, test_case.vload},
                    {I::VTEMP, test_case.vtemp},
                    {I::VLV, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.gastpti,
                                    {{I::VLV, test_case.expected}},
                                    test_case.label);
        }

        // The exhaust-temperature feedback drives the temperature demand.
        setState(fixture.gastpti, {{I::XTEMP, 0.9}, {I::VTEMP, 1.1}});
        success *= (fixture.evaluate() == 0);
        success *= residualsMatch(fixture.gastpti,
                                  {{I::VTEMP, 1.06}},
                                  "temperature feedback");

        // A 1e-4 temperature-gate margin seats the load demand above the
        // temperature demand through the exact inverse smooth ramp, so the
        // smooth LV gate still reproduces the seeded fuel flow.
        auto near_gate_data                   = makeResidualData();
        near_gate_data.parameters[Params::At] = 0.8 + 1.0e-4 / 1.4;
        Fixture<ScalarT> near_gate(near_gate_data);
        near_gate.attachAllInputs();
        success *= near_gate.initialize(0.4);
        success *= stateMatches(near_gate.gastpti,
                                {{I::VLOAD, 0.8155903227031184},
                                 {I::VTEMP, 0.8001000000000001},
                                 {I::VLV, 0.8}},
                                "near-gate initialization");
        success *= scalarMatches(near_gate.input(E::PREF),
                                 0.4077951613515592,
                                 "near-gate published pref");
        success *= (near_gate.evaluate() == 0);
        success *= allResidualsZero(near_gate.gastpti);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// A single rich state and both external inputs drive the two
      /// sensitivity paths; every Enzyme CSR row must match dependency
      /// tracking.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto data = makeResidualData();

        const auto dependency_jacobian = dependencyTrackingJacobian(data, success);
        const auto enzyme_jacobian     = enzymeJacobian(data, success);

        success         *= (dependency_jacobian.size() == enzyme_jacobian.size());
        const auto rows  = std::min(dependency_jacobian.size(), enzyme_jacobian.size());
        for (size_t row = 0; row < rows; ++row)
        {
          if (!isEqual(dependency_jacobian[row], enzyme_jacobian[row], kJacobianTol))
          {
            std::cout << "GASTPTI Jacobian row " << row
                      << " mismatch between dependency tracking and Enzyme\n";
            success = false;
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Governor::GastPtiParameters;
      using Vars   = PhasorDynamics::Governor::GastPtiInternalVariables;
      using Ext    = PhasorDynamics::Governor::GastPtiExternalVariables;
      using Mon    = PhasorDynamics::Governor::GastPtiMonitorableVariables;
      using Mode   = PhasorDynamics::Governor::ResponseMode;
      using Data   = PhasorDynamics::Governor::GastPtiData<RealT, IdxT>;
      using I      = PhasorDynamics::Governor::GastPtiIdx;
      using E      = PhasorDynamics::Governor::GastPtiExt;

      /// A vector row paired with a value: either an input to write or an
      /// expected result. Rows are `GastPtiIdx`/`GastPtiExt` constants, so a
      /// failure report locates itself without any name string to maintain.
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
        std::array<T, E::MAXIMUM>                                   input_values_{};
        std::array<IdxT, E::MAXIMUM>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, E::MAXIMUM> input_nodes_{};

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

          for (size_t port = 0; port < E::MAXIMUM; ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = gastpti.getSignals();
          signals.template attachSignalNode<Ext::OMEGA>(&input_nodes_[E::OMEGA]);
          signals.template attachSignalNode<Ext::PREF>(&input_nodes_[E::PREF]);
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
        fixture.input(E::OMEGA) = static_cast<T>(0.02);
        fixture.input(E::PREF)  = static_cast<T>(0.31);
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Governor::GastPti<T, IdxT>& gastpti) const
      {
        setState(gastpti,
                 {{I::XVALVE, 0.62},
                  {I::XFLOW, 0.55},
                  {I::XTEMP, 0.48},
                  {I::VLOAD, 0.83},
                  {I::VTEMP, 1.35},
                  {I::VLV, 0.71},
                  {I::PMECH, 0.33}});
        setDerivative(gastpti,
                      {{I::XVALVE, 0.011},
                       {I::XFLOW, -0.022},
                       {I::XTEMP, 0.033}});
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

        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.gastpti.y(),
                                   copyVector(explicit_defaults.gastpti.y()),
                                   "documented-default state");
        success *= vectorUnchanged(implicit_defaults.gastpti.yp(),
                                   copyVector(explicit_defaults.gastpti.yp()),
                                   "documented-default derivative");
        success *= vectorUnchanged(implicit_defaults.gastpti.getResidual(),
                                   copyVector(explicit_defaults.gastpti.getResidual()),
                                   "documented-default residual");

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.gastpti);
        setAnswerKeyState(explicit_defaults.gastpti);
        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.gastpti.getResidual(),
                                   copyVector(explicit_defaults.gastpti.getResidual()),
                                   "documented-default dynamic residual");
        return success;
      }

      bool invalidParameterCase(Params parameter, RealT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> model(data);
        return model.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>        unlinked_node;
        PhasorDynamics::Governor::GastPti<ScalarT, IdxT> model(makeData());
        model.getSignals().template attachSignalNode<variable>(&unlinked_node);
        return model.verify() > 0;
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
          success &= rowMatches(static_cast<RealT>(values[i]), snapshot[i], what, i, "changed");
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
                                            const char* label) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        fixture.input(E::PREF) = 77.0; // must stay untouched on rejection
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

        success *= scalarMatches(fixture.pmech(), pmech, "rejected pmech seed preservation");
        success *= scalarMatches(fixture.input(E::OMEGA), 0.0, "rejected omega preservation");
        success *= scalarMatches(fixture.input(E::PREF), 77.0, "rejected pref preservation");
        success *= vectorUnchanged(fixture.gastpti.y(), y_before, "state");
        success *= vectorUnchanged(fixture.gastpti.yp(), yp_before, "derivative");
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
      /// Rows are named by position, which is the `GastPtiIdx` constant the
      /// expectation was written with, leaving no name string to maintain.
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
                  << " mismatch: " << std::setprecision(16) << actual
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
          const auto& [row, expected]  = rows[i];
          success                     &= rowMatches(static_cast<RealT>(values[row]), expected, what, row, context);
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
          success &= rowMatches(static_cast<RealT>(f[row]), 0.0, "residual", row, "at rest");
          success &= rowMatches(static_cast<RealT>(yp[row]), 0.0, "derivative", row, "at rest");
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
        std::cout << label << " mismatch: " << std::setprecision(16) << actual
                  << " != " << expected << "\n";
        return false;
      }

      void noteExpectedLogs(const char* message) const
      {
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << "\n";
        Log::setVerbosity(previous_verbosity);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
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
        for (size_t port = 0; port < E::MAXIMUM; ++port)
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
