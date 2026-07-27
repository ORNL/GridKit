#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
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

      // HYGOV initialization solves the piecewise-linear gate curve exactly
      // while the residual rides its smooth CommonMath approximation. At the
      // mid-segment operating points used here the resulting steady residuals
      // are O(1e-13), so behavioral comparisons use a tolerance well above
      // that gap and well below every pinned answer-key digit.
      static constexpr RealT kBehaviorTol = 1.0e-9;

      // Enzyme and dependency tracking traverse the same smooth expressions
      // differently; their double-precision derivatives agree to O(1e-10).
      static constexpr RealT kJacobianTol = 1.0e-9;

      /// Construction and every verify() error class, including parameter
      /// types, parameter relationships, curve monotonicity, and signal
      /// linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> empty;
        success *= (empty.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (empty.getMonitor() == nullptr);

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> configured(makeData());
        success *= (configured.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() == 0);

        noteExpectedLogs("Testing HYGOV defaults and invalid configurations. "
                         "Logged errors and time-constant warnings are expected.");

        auto minimal_data                      = makeMinimalData();
        minimal_data.parameters[Params::Trate] = 100.0;
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> minimal(minimal_data);
        success *= (minimal.verify() == 0);
        success *= defaultsMatchDocumentedValues();

        success *= (empty.verify() > 0);

        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> missing_trate(makeMinimalData());
        success *= (missing_trate.verify() > 0);

        success *= invalidParameterCase(Params::Trate, 0.0);
        success *= invalidParameterCase(Params::Rtemp, 0.0);
        success *= invalidParameterCase(Params::Tr, -0.1);
        success *= invalidParameterCase(Params::Tf, -0.1);
        success *= invalidParameterCase(Params::Tg, -0.1);
        success *= invalidParameterCase(Params::Tw, -0.1);
        success *= invalidParameterCase(Params::Tn, -0.1);
        success *= invalidParameterCase(Params::Tnp, -0.1);
        success *= invalidParameterCase(Params::Velm, -0.1);
        success *= invalidParameterCase(Params::Gmin, 1.1);
        success *= invalidParameterCase(Params::At, 0.0);
        success *= invalidParameterCase(Params::Dturb, -0.1);
        success *= invalidParameterCase(Params::db1, -0.1);
        success *= invalidParameterCase(Params::Hdam, 0.0);
        success *= invalidParameterCase(Params::Gv2, 0.1);
        success *= invalidParameterCase(Params::Pgv2, 0.1);

        // db2 is accepted for source-format compatibility and never used.
        auto backlash_data                    = makeData();
        backlash_data.parameters[Params::db2] = 0.5;
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> backlash_model(backlash_data);
        success *= (backlash_model.verify() == 0);

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                   = makeData();
        integer_real.parameters[Params::Tw] = static_cast<IdxT>(2);
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> integer_real_model(integer_real);
        success *= (integer_real_model.verify() == 0);

        auto bad_numeric_type                      = makeData();
        bad_numeric_type.parameters[Params::Trate] = true;
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> bad_numeric_model(bad_numeric_type);
        success *= (bad_numeric_model.verify() > 0);

        success *= unlinkedSignalRejected<Ext::OMEGA>();
        success *= unlinkedSignalRejected<Ext::PREF>();
        success *= unlinkedSignalRejected<Ext::PAUX>();

        // All five zero time constants use the documented numerical floor and
        // still admit a consistent steady-state initialization.
        auto zero_time                    = makeData();
        zero_time.parameters[Params::Tr]  = 0.0;
        zero_time.parameters[Params::Tf]  = 0.0;
        zero_time.parameters[Params::Tg]  = 0.0;
        zero_time.parameters[Params::Tw]  = 0.0;
        zero_time.parameters[Params::Tnp] = 0.0;

        Fixture<ScalarT> fixture(zero_time);
        success *= fixture.initialize(0.4);
        success *= (fixture.evaluate() == 0);
        success *= allResidualsZero(fixture.hygov);

        return success.report(__func__);
      }

      /// A nonidentity power-base initialization with every port attached.
      /// The machine-seeded pmech node must remain unchanged while HYGOV
      /// initializes and publishes its resolved load reference.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        auto data                      = makeData();
        data.parameters[Params::Trate] = 50.0;

        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        fixture.input(E::PAUX)  = 0.02;
        fixture.input(E::PREF)  = 99.0; // stale value the publication must replace
        success                *= fixture.initialize(0.4);
        success                *= (fixture.hygov.tagDifferentiable() == 0);
        success                *= (fixture.evaluate() == 0);

        const auto* y  = fixture.hygov.y().getData();
        success       *= scalarMatches(y[I::XF], 0.0, "XF at rest");
        success       *= scalarMatches(y[I::C], 0.9, "C on component base");
        success       *= scalarMatches(y[I::G], 0.9, "G on component base");
        success       *= scalarMatches(y[I::Q], 0.9, "Q on component base");
        success       *= scalarMatches(y[I::PGV], 0.9, "PGV on component base");
        success       *= scalarMatches(y[I::H], 1.0, "H at the dam head");
        success       *= scalarMatches(fixture.pmech(), 0.4, "preserved pmech seed");

        success *= scalarMatches(fixture.input(E::OMEGA), 0.0, "preserved omega input");
        success *= scalarMatches(fixture.input(E::PREF), 0.0025, "published pref");
        success *= scalarMatches(fixture.input(E::PAUX), 0.02, "preserved paux input");

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
          success *= scalarMatches(monitored[3], 0.9, "monitored desiredgate");
          success *= scalarMatches(monitored[4], 0.9, "monitored gate");
          success *= scalarMatches(monitored[5], 0.9, "monitored flow");
          success *= scalarMatches(monitored[6], 1.0, "monitored head");
        }
        else
        {
          std::cout << "HYGOV monitor emitted " << monitored.size()
                    << " values instead of 7\n";
          success = false;
        }

        for (size_t i = 0; i < static_cast<size_t>(fixture.hygov.size()); ++i)
        {
          const bool expected = i <= I::Q;
          if (fixture.hygov.tag()[i] != expected)
          {
            std::cout << "HYGOV differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= allResidualsZero(fixture.hygov);

        // A system-base reference step lands on the governor error scaled by
        // the base ratio.
        fixture.input(E::PREF)  = 0.1025; // the published 0.0025 plus a 0.1 step
        success                *= (fixture.evaluate() == 0);
        success                *= residualsMatch(fixture.hygov,
                                                 {{I::EF, 0.2}},
                                  "reference step on the component base");

        // Unattached ports fall back to the references latched by
        // initialize(), so the same steady state holds without a controller.
        Fixture<ScalarT> fallback(data);
        success *= fallback.initialize(0.4);
        success *= (fallback.evaluate() == 0);
        success *= allResidualsZero(fallback.hygov);

        return success.report(__func__);
      }

      /// Mechanical-power and gate-limit initialization domains. Every
      /// rejection is atomic; exact Gmin/Gmax boundaries and a zero power
      /// seed remain admissible.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible HYGOV mechanical-power and gate "
                         "initialization points. Logged errors are expected.");

        struct RejectionCase
        {
          const char* label;
          RealT       pmech;
          RealT       gmin;
          RealT       gmax;
        };

        const std::array<RejectionCase, 4> rejected{{
            {"mechanical power above the gate curve", 1.0, 0.05, 0.95},
            {"mechanical power below the gate curve", -0.3, 0.05, 0.95},
            {"initialized gate above Gmax", 0.4, 0.05, 0.5},
            {"initialized gate below Gmin", 0.4, 0.6, 0.95},
        }};

        for (const auto& test_case : rejected)
        {
          auto data                      = makeResidualData();
          data.parameters[Params::Gmin]  = test_case.gmin;
          data.parameters[Params::Gmax]  = test_case.gmax;
          success                       *= initializationRejectedAtomically(
              data, test_case.pmech, test_case.label);
        }

        // An invalid configuration is rejected before any state is written.
        auto invalid_data                      = makeResidualData();
        invalid_data.parameters[Params::Rtemp] = 0.0;
        Fixture<ScalarT> invalid_fixture(invalid_data);
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

        // Exact gate-limit boundaries and a zero power seed stay admissible.
        struct AdmissibleCase
        {
          const char* label;
          RealT       pmech;
          RealT       gmin;
          RealT       gmax;
          RealT       gate;
        };

        for (const auto& accepted : std::array<AdmissibleCase, 3>{{
                 {"gate landing exactly on Gmax", 0.8, 0.0, 0.9, 0.9},
                 {"gate landing exactly on Gmin", 0.0, 0.1, 1.0, 0.1},
                 {"zero mechanical-power seed", 0.0, 0.0, 1.0, 0.1},
             }})
        {
          auto data                     = makeData();
          data.parameters[Params::Gmin] = accepted.gmin;
          data.parameters[Params::Gmax] = accepted.gmax;

          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize(accepted.pmech);
          success *= stateMatches(fixture.hygov,
                                  {{I::C, accepted.gate}, {I::G, accepted.gate}},
                                  accepted.label);
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

        // Values are pinned after an independent one-time evaluation of the
        // documented equations at setAnswerKeyState()/setAnswerKeyInputs().
        const std::array<Row, I::MAXIMUM> expected{{
            {I::XN, -0.07785714285714286},
            {I::XF, -0.7300000000000001},
            {I::C, 0.06},
            {I::G, 0.1233333333333334},
            {I::Q, 0.011538461538461414},
            {I::OMEGADB, 0.0033514666467982894},
            {I::EF, 0.5863},
            {I::FC, -1.8512500000000003},
            {I::RC, 0.029996890386450745},
            {I::PGV, -0.04600000003160343},
            {I::H, -0.033299999999999885},
            {I::PMECH, -0.012679999999999934},
        }};

        success *= (static_cast<size_t>(fixture.hygov.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.hygov, expected);

        return success.report(__func__);
      }

      /// Speed deadband, desired-gate velocity limiting, gate-position
      /// anti-windup, and turbine damping at nonzero speed deviation.
      TestOutcome governorControl()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);

        // The type-1 deadband below, inside, and above the +-0.01 band.
        struct DeadbandCase
        {
          RealT omega;
          RealT expected;
        };

        for (const auto& test_case : std::array<DeadbandCase, 3>{{
                 {-0.05, -0.049996641662021946},
                 {0.004, 0.0009004582873718001},
                 {0.05, 0.049996641662021946},
             }})
        {
          fixture.input(E::OMEGA) = test_case.omega;
          setState(fixture.hygov, {{I::OMEGADB, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.hygov,
                                    {{I::OMEGADB, test_case.expected}},
                                    "speed deadband");
        }
        fixture.input(E::OMEGA) = 0.0;

        // The desired-gate velocity target driven below, inside, and above
        // the +-Velm rate limit.
        struct VelocityCase
        {
          RealT fc;
          RealT expected;
        };

        for (const auto& test_case : std::array<VelocityCase, 3>{{
                 {-0.6, -0.15},
                 {0.05, 0.04999999999984272},
                 {0.6, 0.15000000000000002},
             }})
        {
          setState(fixture.hygov, {{I::FC, test_case.fc}, {I::RC, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.hygov,
                                    {{I::RC, test_case.expected}},
                                    "gate velocity limit");
        }

        // The desired-gate anti-windup at three controller directions: both
        // saturations block an outward rate and Gmax admits a restoring one.
        struct AntiWindupCase
        {
          const char* label;
          RealT       c;
          RealT       rc;
          RealT       expected;
        };

        for (const auto& test_case : std::array<AntiWindupCase, 3>{{
                 {"Gmax blocks an outward desired-gate rate", 1.2, 0.2, 0.0},
                 {"Gmin blocks an outward desired-gate rate", -0.2, -0.2, 0.0},
                 {"Gmax admits a restoring desired-gate rate", 1.2, -0.2, -0.2},
             }})
        {
          setState(fixture.hygov, {{I::C, test_case.c}, {I::RC, test_case.rc}});
          setDerivative(fixture.hygov, {{I::C, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.hygov,
                                    {{I::C, test_case.expected}},
                                    test_case.label);
        }

        // Turbine damping proportional to speed deviation and gate.
        fixture.input(E::OMEGA) = 0.05;
        setState(fixture.hygov,
                 {{I::G, 0.6}, {I::Q, 0.7}, {I::H, 1.1}, {I::PMECH, 0.5}});
        success *= (fixture.evaluate() == 0);
        success *= residualsMatch(fixture.hygov,
                                  {{I::PMECH, -0.2677999999999999}},
                                  "turbine damping");

        return success.report(__func__);
      }

      /// The nonlinear gate-power curve on every rising segment, the water
      /// column away from the dam head, and initialization through the curve
      /// with and without the speed-damping term.
      TestOutcome turbineDynamics()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData());
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);

        // One gate point inside each of the five curve segments.
        struct CurveCase
        {
          RealT gate;
          RealT expected;
        };

        for (const auto& test_case : std::array<CurveCase, 5>{{
                 {0.1, 0.07500000000021236},
                 {0.3, 0.28500000000007075},
                 {0.5, 0.5399999999999371},
                 {0.7, 0.7549999999999292},
                 {0.9, 0.9249999999998506},
             }})
        {
          setState(fixture.hygov, {{I::G, test_case.gate}, {I::PGV, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.hygov,
                                    {{I::PGV, test_case.expected}},
                                    "gate-power curve");
        }

        // A head away from the dam head drives the flow and head rows.
        setState(fixture.hygov, {{I::Q, 0.61}, {I::H, 0.9}, {I::PGV, 0.55}});
        setDerivative(fixture.hygov, {{I::Q, 0.05}});
        success *= (fixture.evaluate() == 0);
        success *= residualsMatch(fixture.hygov,
                                  {{I::Q, 0.18076923076923068}, {I::H, -0.09984999999999994}},
                                  "water column");

        // A seed inside the third curve segment initializes through the
        // nonidentity curve inversion.
        Fixture<ScalarT> curve_fixture(makeResidualData());
        curve_fixture.attachAllInputs();
        success *= curve_fixture.initialize(0.33761676);
        success *= stateMatches(curve_fixture.hygov,
                                {{I::C, 0.5000001394782843}, {I::G, 0.5000001394782843}},
                                "nonidentity curve inversion");
        success *= scalarMatches(curve_fixture.input(E::PREF),
                                 0.015000004184348527,
                                 "nonidentity-curve published pref");
        success *= scalarMatches(curve_fixture.pmech(), 0.33761676, "preserved pmech seed");
        success *= (curve_fixture.evaluate() == 0);
        success *= allResidualsZero(curve_fixture.hygov);

        // A nonzero initial speed deviation folds the Dturb damping loss into
        // the gate solve, so the damped point still initializes at rest.
        Fixture<ScalarT> damped_fixture(makeResidualData());
        damped_fixture.attachAllInputs();
        damped_fixture.input(E::OMEGA)  = 0.03;
        success                        *= damped_fixture.initialize(0.48676047);
        success                        *= stateMatches(damped_fixture.hygov,
                                                       {{I::G, 0.7000002496006894},
                                                        {I::OMEGADB, 0.029757154589893794}},
                                "damped initialization");
        success                        *= scalarMatches(damped_fixture.input(E::PREF),
                                 0.03587858478296758,
                                 "damped published pref");
        success                        *= (damped_fixture.evaluate() == 0);
        success                        *= allResidualsZero(damped_fixture.hygov);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// A single rich state and all three external inputs drive both
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
            std::cout << "HYGOV Jacobian row " << row
                      << " mismatch between dependency tracking and Enzyme\n";
            success = false;
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Governor::HygovParameters;
      using Vars   = PhasorDynamics::Governor::HygovInternalVariables;
      using Ext    = PhasorDynamics::Governor::HygovExternalVariables;
      using Mon    = PhasorDynamics::Governor::HygovMonitorableVariables;
      using Data   = PhasorDynamics::Governor::HygovData<RealT, IdxT>;
      using I      = PhasorDynamics::Governor::HygovIdx;
      using E      = PhasorDynamics::Governor::HygovExt;

      /// A vector row paired with a value: either an input to write or an
      /// expected result. Rows are `HygovIdx`/`HygovExt` constants, so a
      /// failure report locates itself without any name string to maintain.
      using Row    = std::pair<size_t, RealT>;
      using Rows   = std::initializer_list<Row>;
      using HygovT = PhasorDynamics::Governor::Hygov<ScalarT, IdxT>;

      /// Owns the HYGOV model, the assigned mechanical-power node, and the
      /// attached input nodes. Signal storage is declared before the model so
      /// every referenced node outlives HYGOV. Copying would invalidate the
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
          : hygov(data)
        {
          hygov.setSystemBase(60.0, system_va_base);
          hygov.getSignals().template assignSignalNode<Vars::PMECH>(&pmech_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Attach fixture-owned storage to every external input.
        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = hygov.size();

          for (size_t port = 0; port < E::MAXIMUM; ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = hygov.getSignals();
          signals.template attachSignalNode<Ext::OMEGA>(&input_nodes_[E::OMEGA]);
          signals.template attachSignalNode<Ext::PREF>(&input_nodes_[E::PREF]);
          signals.template attachSignalNode<Ext::PAUX>(&input_nodes_[E::PAUX]);
        }

        /// Seed the assigned mechanical-power node on the system base.
        void seedPmech(RealT pmech)
        {
          pmech_node_.init(static_cast<T>(pmech));
        }

        /// Everything HYGOV initialization requires: allocation,
        /// verification, and a machine-seeded mechanical-power node.
        bool prepare(RealT pmech)
        {
          const bool success = (hygov.allocate() == 0) && (hygov.verify() == 0);
          if (!success)
          {
            std::cout << "HYGOV fixture preparation failed\n";
            return false;
          }

          seedPmech(pmech);
          return true;
        }

        /// prepare() plus successful HYGOV initialization.
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

        T& input(size_t port)
        {
          return input_values_[port];
        }

        IdxT inputIndex(size_t port) const
        {
          return input_indices_[port];
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
        auto data = makeMinimalData();

        // These are the documented defaults. The all-zero source curve
        // selects the identity curve, spelled out here point by point.
        data.parameters[Params::Trate] = 100.0;
        data.parameters[Params::Rperm] = 0.04;
        data.parameters[Params::Rtemp] = 0.3;
        data.parameters[Params::Tr]    = 5.0;
        data.parameters[Params::Tf]    = 0.05;
        data.parameters[Params::Tg]    = 0.5;
        data.parameters[Params::Velm]  = 0.2;
        data.parameters[Params::Gmax]  = 1.0;
        data.parameters[Params::Gmin]  = 0.0;
        data.parameters[Params::Tw]    = 1.0;
        data.parameters[Params::At]    = 1.2;
        data.parameters[Params::Dturb] = 0.5;
        data.parameters[Params::Qnl]   = 0.05;
        data.parameters[Params::Tn]    = 0.0;
        data.parameters[Params::Tnp]   = 0.0;
        data.parameters[Params::db1]   = 0.0;
        data.parameters[Params::db2]   = 0.0;
        data.parameters[Params::Hdam]  = 1.0;
        data.parameters[Params::Gv0]   = 0.0;
        data.parameters[Params::Gv1]   = 0.2;
        data.parameters[Params::Gv2]   = 0.4;
        data.parameters[Params::Gv3]   = 0.6;
        data.parameters[Params::Gv4]   = 0.8;
        data.parameters[Params::Gv5]   = 1.0;
        data.parameters[Params::Pgv0]  = 0.0;
        data.parameters[Params::Pgv1]  = 0.2;
        data.parameters[Params::Pgv2]  = 0.4;
        data.parameters[Params::Pgv3]  = 0.6;
        data.parameters[Params::Pgv4]  = 0.8;
        data.parameters[Params::Pgv5]  = 1.0;
        return data;
      }

      Data makeData() const
      {
        auto data = makeMinimalData();

        data.parameters[Params::Trate] = 100.0;
        data.parameters[Params::Rperm] = 0.05;
        data.parameters[Params::Rtemp] = 0.4;
        data.parameters[Params::Tr]    = 5.0;
        data.parameters[Params::Tf]    = 0.2;
        data.parameters[Params::Tg]    = 0.5;
        data.parameters[Params::Velm]  = 0.5;
        data.parameters[Params::Gmax]  = 1.0;
        data.parameters[Params::Gmin]  = 0.0;
        data.parameters[Params::Tw]    = 1.0;
        data.parameters[Params::At]    = 1.0;
        data.parameters[Params::Dturb] = 0.0;
        data.parameters[Params::Qnl]   = 0.1;
        data.parameters[Params::Tn]    = 0.0;
        data.parameters[Params::Tnp]   = 1.0;
        data.parameters[Params::db1]   = 0.0;
        data.parameters[Params::db2]   = 0.0;
        data.parameters[Params::Hdam]  = 1.0;
        data.parameters[Params::Gv0]   = 0.0;
        data.parameters[Params::Gv1]   = 0.2;
        data.parameters[Params::Gv2]   = 0.4;
        data.parameters[Params::Gv3]   = 0.6;
        data.parameters[Params::Gv4]   = 0.8;
        data.parameters[Params::Gv5]   = 1.0;
        data.parameters[Params::Pgv0]  = 0.0;
        data.parameters[Params::Pgv1]  = 0.2;
        data.parameters[Params::Pgv2]  = 0.4;
        data.parameters[Params::Pgv3]  = 0.6;
        data.parameters[Params::Pgv4]  = 0.8;
        data.parameters[Params::Pgv5]  = 1.0;
        return data;
      }

      Data makeResidualData() const
      {
        auto data = makeData();

        data.parameters[Params::Trate] = 50.0;
        data.parameters[Params::Rperm] = 0.06;
        data.parameters[Params::Rtemp] = 0.4;
        data.parameters[Params::Tr]    = 4.0;
        data.parameters[Params::Tf]    = 0.2;
        data.parameters[Params::Tg]    = 0.6;
        data.parameters[Params::Velm]  = 0.15;
        data.parameters[Params::Gmax]  = 0.95;
        data.parameters[Params::Gmin]  = 0.05;
        data.parameters[Params::Tw]    = 1.3;
        data.parameters[Params::At]    = 1.1;
        data.parameters[Params::Dturb] = 0.6;
        data.parameters[Params::Qnl]   = 0.08;
        data.parameters[Params::Tn]    = 0.7;
        data.parameters[Params::Tnp]   = 1.4;
        data.parameters[Params::db1]   = 0.01;
        data.parameters[Params::Hdam]  = 1.2;
        data.parameters[Params::Pgv1]  = 0.15;
        data.parameters[Params::Pgv2]  = 0.42;
        data.parameters[Params::Pgv3]  = 0.66;
        data.parameters[Params::Pgv4]  = 0.85;
        return data;
      }

      /// The external inputs the residual answer key is evaluated against.
      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(E::OMEGA) = static_cast<T>(0.02);
        fixture.input(E::PREF)  = static_cast<T>(0.31);
        fixture.input(E::PAUX)  = static_cast<T>(0.07);
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Governor::Hygov<T, IdxT>& hygov) const
      {
        setState(hygov,
                 {{I::XN, 0.11},
                  {I::XF, 0.23},
                  {I::C, 0.52},
                  {I::G, 0.47},
                  {I::Q, 0.61},
                  {I::OMEGADB, 0.015},
                  {I::EF, 0.08},
                  {I::FC, 0.12},
                  {I::RC, 0.09},
                  {I::PGV, 0.55},
                  {I::H, 1.12},
                  {I::PMECH, 0.33}});
        setDerivative(hygov,
                      {{I::XN, 0.01},
                       {I::XF, -0.02},
                       {I::C, 0.03},
                       {I::G, -0.04},
                       {I::Q, 0.05}});
      }

      /// Omitting every optional parameter must give exactly the model built
      /// from the defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        auto implicit_data                      = makeMinimalData();
        implicit_data.parameters[Params::Trate] = 100.0;

        Fixture<ScalarT> implicit_defaults(implicit_data);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData());
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        bool success = implicit_defaults.initialize(0.3)
                       && explicit_defaults.initialize(0.3);
        if (!success)
        {
          std::cout << "HYGOV documented-default comparison failed to initialize\n";
          return false;
        }

        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.hygov.y(),
                                   copyVector(explicit_defaults.hygov.y()),
                                   "documented-default state");
        success *= vectorUnchanged(implicit_defaults.hygov.yp(),
                                   copyVector(explicit_defaults.hygov.yp()),
                                   "documented-default derivative");
        success *= vectorUnchanged(implicit_defaults.hygov.getResidual(),
                                   copyVector(explicit_defaults.hygov.getResidual()),
                                   "documented-default residual");

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.hygov);
        setAnswerKeyState(explicit_defaults.hygov);
        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.hygov.getResidual(),
                                   copyVector(explicit_defaults.hygov.getResidual()),
                                   "documented-default dynamic residual");
        return success;
      }

      bool invalidParameterCase(Params parameter, RealT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> model(data);
        return model.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>      unlinked_node;
        PhasorDynamics::Governor::Hygov<ScalarT, IdxT> model(makeData());
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
        auto* y  = fixture.hygov.y().getData();
        auto* yp = fixture.hygov.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.hygov.y().getSize()); ++i)
        {
          y[i]  = 0.125 + 0.01 * static_cast<RealT>(i);
          yp[i] = -0.25 - 0.01 * static_cast<RealT>(i);
        }
        fixture.seedPmech(pmech);
        fixture.hygov.y().setDataUpdated();
        fixture.hygov.yp().setDataUpdated();
      }

      bool initializationRejectedAtomically(const Data& data,
                                            RealT       pmech,
                                            const char* label) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        fixture.input(E::PAUX) = 0.02;
        fixture.input(E::PREF) = 77.0; // must stay untouched on rejection
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

        success *= scalarMatches(fixture.pmech(), pmech, "rejected pmech seed preservation");
        success *= scalarMatches(fixture.input(E::OMEGA), 0.0, "rejected omega preservation");
        success *= scalarMatches(fixture.input(E::PREF), 77.0, "rejected pref preservation");
        success *= scalarMatches(fixture.input(E::PAUX), 0.02, "rejected paux preservation");
        success *= vectorUnchanged(fixture.hygov.y(), y_before, "state");
        success *= vectorUnchanged(fixture.hygov.yp(), yp_before, "derivative");
        return success;
      }

      /// Write state rows and publish the update, folding in the
      /// setDataUpdated() that a hand-written write block has to remember.
      template <typename T>
      void setState(PhasorDynamics::Governor::Hygov<T, IdxT>& hygov, Rows rows) const
      {
        auto* y = hygov.y().getData();
        for (const auto& [row, value] : rows)
        {
          y[row] = static_cast<T>(value);
        }
        hygov.y().setDataUpdated();
      }

      /// setState() for the derivative vector.
      template <typename T>
      void setDerivative(PhasorDynamics::Governor::Hygov<T, IdxT>& hygov, Rows rows) const
      {
        auto* yp = hygov.yp().getData();
        for (const auto& [row, value] : rows)
        {
          yp[row] = static_cast<T>(value);
        }
        hygov.yp().setDataUpdated();
      }

      /// Compare one vector row against its expected value. Every row check
      /// in this suite reports through here, so failures share one format.
      /// Rows are named by position, which is the `HygovIdx` constant the
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
        std::cout << "HYGOV " << what << " row " << row << ' ' << context
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

      bool residualsMatch(const HygovT& hygov, Rows rows, const char* context = "") const
      {
        return rowsMatch(hygov.getResidual(), rows.begin(), rows.size(), "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const HygovT&                hygov,
                          const std::array<Row, size>& rows,
                          const char*                  context = "") const
      {
        return rowsMatch(hygov.getResidual(), rows.data(), size, "residual", context);
      }

      bool stateMatches(const HygovT& hygov, Rows rows, const char* context = "") const
      {
        return rowsMatch(hygov.y(), rows.begin(), rows.size(), "state", context);
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
        auto* y  = fixture.hygov.y().getData();
        auto* yp = fixture.hygov.yp().getData();

        const auto model_size = static_cast<size_t>(fixture.hygov.size());
        for (size_t i = 0; i < model_size; ++i)
        {
          y[i].setVariableNumber(i);
          yp[i].setVariableNumber(i);
        }
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          fixture.input(port).setVariableNumber(fixture.inputIndex(port));
        }

        fixture.hygov.y().setDataUpdated();
        fixture.hygov.yp().setDataUpdated();
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
        setAnswerKeyState(fixture.hygov);
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
          TestStatus& success) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.4);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.hygov);
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
