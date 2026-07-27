#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
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
    class ExciterEsdc1aTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ExciterEsdc1aTests()  = default;
      ~ExciterEsdc1aTests() = default;

      // ESDC1A initialization seeds the smooth high-value gate through the
      // ramp inverse, leaving steady residuals of O(1e-12). Behavioral
      // comparisons use a tolerance three orders above that guard.
      static constexpr RealT kBehaviorTol = 1.0e-9;

      // Enzyme and dependency tracking traverse the same smooth expressions
      // differently; their double-precision derivatives agree to O(1e-10).
      static constexpr RealT kJacobianTol = 1.0e-9;

      /// Construction and every verify() error class, including parameter
      /// types, parameter relationships, bus ownership, and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing ESDC1A defaults and invalid configurations. "
                         "Logged errors and time-constant warnings are expected.");

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Exciter::Esdc1a<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (empty.getMonitor() == nullptr);
        success *= (empty.verify() > 0);

        Fixture<ScalarT> configured(makeData());
        success *= (configured.esdc1a.size() == static_cast<IdxT>(I::MAXIMUM));
        success *= (configured.esdc1a.getMonitor() != nullptr);
        success *= configured.prepare(1.2);

        // ESDC1A has no required parameters: an empty parameter set is the
        // documented-default model.
        Fixture<ScalarT> minimal(makeMinimalData());
        success *= minimal.prepare(1.2);
        success *= defaultsMatchDocumentedValues();

        // A model without the required efd output assignment is rejected.
        PhasorDynamics::Exciter::Esdc1a<ScalarT, IdxT> unassigned(&bus, makeData());
        success *= (unassigned.verify() > 0);

        success *= invalidParameterCase(Params::Ka, 0.0);
        success *= invalidParameterCase(Params::Ta, 0.0);
        success *= invalidParameterCase(Params::Te, 0.0);
        success *= invalidParameterCase(Params::Tc, -0.1);
        success *= invalidParameterCase(Params::Tr, -0.1);
        success *= invalidParameterCase(Params::Tb, -0.1);
        success *= invalidParameterCase(Params::Tf1, -0.1);
        success *= invalidParameterCase(Params::Vrmin, 2.0);
        success *= invalidParameterCase(Params::UEL, static_cast<IdxT>(4));
        success *= invalidParameterCase(Params::UEL, static_cast<RealT>(2.5));
        success *= invalidParameterCase(Params::UEL, true);
        success *= invalidParameterCase(Params::Se1, 0.0);
        success *= invalidParameterCase(Params::E2, 2.8);
        success *= invalidParameterCase(Params::Se2, 0.08);
        success *= invalidParameterCase(Params::E1, -1.0);

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                   = makeData();
        integer_real.parameters[Params::Ka] = static_cast<IdxT>(40);
        Fixture<ScalarT> integer_real_fixture(integer_real);
        success *= (integer_real_fixture.esdc1a.verify() == 0);
        success *= invalidParameterCase(Params::Ka, true);

        for (const Params flag : {Params::Spdmlt, Params::exclim})
        {
          auto bad_integer             = makeData();
          bad_integer.parameters[flag] = static_cast<IdxT>(2);
          Fixture<ScalarT> bad_integer_fixture(bad_integer);
          success *= (bad_integer_fixture.esdc1a.verify() > 0);

          auto bad_real             = makeData();
          bad_real.parameters[flag] = static_cast<RealT>(0.5);
          Fixture<ScalarT> bad_real_fixture(bad_real);
          success *= (bad_real_fixture.esdc1a.verify() > 0);
        }

        // Real-valued 0/1, integer 0/1, and JSON booleans are all accepted
        // for the two switches, matching ESDC1A's REPCA-style parameter
        // contract.
        auto switch_forms                       = makeData();
        switch_forms.parameters[Params::Spdmlt] = static_cast<RealT>(1.0);
        switch_forms.parameters[Params::exclim] = static_cast<IdxT>(0);
        Fixture<ScalarT> switch_fixture(switch_forms);
        switch_fixture.attachAllInputs();
        success *= (switch_fixture.esdc1a.verify() == 0);

        // The enabled speed multiplier requires an attached speed input.
        auto speed_required                       = makeData();
        speed_required.parameters[Params::Spdmlt] = true;
        Fixture<ScalarT> speed_required_fixture(speed_required);
        success *= (speed_required_fixture.esdc1a.verify() > 0);

        PhasorDynamics::SignalNode<ScalarT, IdxT>      busless_efd_node;
        PhasorDynamics::Exciter::Esdc1a<ScalarT, IdxT> busless(nullptr, makeData());
        busless.getSignals().template assignSignalNode<Vars::EFD>(&busless_efd_node);
        success *= (busless.verify() > 0);

        success *= unlinkedSignalRejected<Ext::OMEGA>();
        success *= unlinkedSignalRejected<Ext::VREF>();
        success *= unlinkedSignalRejected<Ext::VS>();
        success *= unlinkedSignalRejected<Ext::VUEL>();

        // All three floored time constants at zero use the documented
        // numerical floor and still admit a consistent steady-state
        // initialization.
        auto zero_time                    = makeData();
        zero_time.parameters[Params::Tr]  = 0.0;
        zero_time.parameters[Params::Tb]  = 0.0;
        zero_time.parameters[Params::Tf1] = 0.0;

        Fixture<ScalarT> floored(zero_time);
        success *= floored.initialize(1.2);
        success *= (floored.evaluate() == 0);
        success *= allResidualsZero(floored.esdc1a);

        return success.report(__func__);
      }

      /// Initialization with every port attached: the seeded field voltage
      /// is preserved, the resolved voltage reference is published, the
      /// known inputs are read but never overwritten, and every selector
      /// combination reaches a zero-derivative, zero-residual state.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeData());
        fixture.attachAllInputs(99.0);
        fixture.input(E::OMEGA)  = 0.02;
        fixture.input(E::VS)     = 0.03;
        fixture.input(E::VUEL)   = -0.4;
        success                 *= fixture.initialize(1.2);
        success                 *= (fixture.esdc1a.tagDifferentiable() == 0);
        success                 *= (fixture.evaluate() == 0);

        const auto* y  = fixture.esdc1a.y().getData();
        success       *= scalarMatches(y[I::EFDP], 1.2, "EFDP");
        success       *= scalarMatches(y[I::VC], 1.0, "VC");
        success       *= scalarMatches(y[I::VR], 0.12, "VR");
        success       *= scalarMatches(y[I::VF], 0.0, "VF");
        success       *= scalarMatches(y[I::EV], 0.003, "EV gate input");
        success       *= scalarMatches(y[I::VHV], 0.003, "VHV");
        success       *= scalarMatches(y[I::SE], 0.0, "SE");
        success       *= scalarMatches(y[I::VFE], 0.12, "VFE");
        success       *= scalarMatches(fixture.efd(), 1.2, "seeded efd");

        success *= scalarMatches(fixture.input(E::VREF), 0.973, "published vref");
        success *= scalarMatches(fixture.input(E::OMEGA), 0.02, "preserved speed input");
        success *= scalarMatches(fixture.input(E::VS), 0.03, "preserved vs input");
        success *= scalarMatches(fixture.input(E::VUEL), -0.4, "preserved vuel input");

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.esdc1a.getMonitor());
        std::stringstream monitor_output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, monitor_output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string monitor_header;
        std::string monitor_values;
        std::getline(monitor_output, monitor_header);
        std::getline(monitor_output, monitor_values);
        success              *= (monitor_header == "t,Esdc1a_esdc1a_test_efd,Esdc1a_esdc1a_test_vc,"
                                                   "Esdc1a_esdc1a_test_vr,Esdc1a_esdc1a_test_vf,"
                                                   "Esdc1a_esdc1a_test_se,Esdc1a_esdc1a_test_vfe");
        const auto monitored  = Tokenizer<RealT>(monitor_values, ',')();
        if (monitored.size() == 7)
        {
          success *= scalarMatches(monitored[1], 1.2, "monitored efd");
          success *= scalarMatches(monitored[2], 1.0, "monitored vc");
          success *= scalarMatches(monitored[3], 0.12, "monitored vr");
          success *= scalarMatches(monitored[4], 0.0, "monitored vf");
          success *= scalarMatches(monitored[5], 0.0, "monitored se");
          success *= scalarMatches(monitored[6], 0.12, "monitored vfe");
        }
        else
        {
          std::cout << "ESDC1A monitor emitted " << monitored.size()
                    << " values instead of 7\n";
          success = false;
        }

        for (size_t i = 0; i < static_cast<size_t>(fixture.esdc1a.size()); ++i)
        {
          const bool expected = i <= I::XLL;
          if (fixture.esdc1a.tag()[i] != expected)
          {
            std::cout << "ESDC1A differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= allResidualsZero(fixture.esdc1a);

        // With no attached inputs the latched values act as constant
        // references.
        Fixture<ScalarT> unattached(makeData());
        success *= unattached.initialize(1.2);
        success *= (unattached.evaluate() == 0);
        success *= allResidualsZero(unattached.esdc1a);

        // Every selector combination must preserve the seeded field voltage
        // and produce a zero-derivative, zero-residual state.
        for (IdxT uel = 0; uel < 4; ++uel)
        {
          for (const bool speed_flag : {false, true})
          {
            for (const bool limit_flag : {false, true})
            {
              auto scenario_data                       = makeData();
              scenario_data.parameters[Params::UEL]    = uel;
              scenario_data.parameters[Params::Spdmlt] = speed_flag;
              scenario_data.parameters[Params::exclim] = limit_flag;

              Fixture<ScalarT> scenario(scenario_data);
              scenario.attachAllInputs();
              scenario.input(E::OMEGA) = 0.02;
              scenario.input(E::VS)    = 0.03;
              scenario.input(E::VUEL)  = -0.4;
              if (!scenario.initialize(1.2))
              {
                std::cout << "ESDC1A initialization scenario failed: UEL=" << uel
                          << ", Spdmlt=" << speed_flag
                          << ", exclim=" << limit_flag << '\n';
                success = false;
                continue;
              }

              success *= (scenario.evaluate() == 0);
              success *= allResidualsZero(scenario.esdc1a);
              success *= scalarMatches(scenario.efd(), 1.2, "scenario efd preservation");
            }
          }
        }

        return success.report(__func__);
      }

      /// The inadmissible initialization points: every rejection is atomic,
      /// and the admissible operating points next to them still initialize
      /// to zero residuals.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible ESDC1A initialization points. "
                         "Logged errors are expected.");

        // An enabled speed multiplier with omega = -1 zeroes the seed
        // denominator.
        auto speed_data                        = makeData();
        speed_data.parameters[Params::Spdmlt]  = true;
        success                               *= initializationRejectedAtomically(speed_data,
                                                    1.2,
                                                                                  {-1.0, 77.0, 77.0, -77.0},
                                                    "zero speed-multiplier denominator");

        // The seeded field voltage maps to a regulator output above Vrmax.
        auto limit_data                       = makeData();
        limit_data.parameters[Params::Vrmax]  = 0.05;
        limit_data.parameters[Params::Vrmin]  = -0.05;
        success                              *= initializationRejectedAtomically(limit_data,
                                                    1.2,
                                                                                 {0.0, 77.0, 77.0, -77.0},
                                                    "regulator output outside limits");

        // A UEL input above the gate operating point holds the high-value
        // gate active, which the smooth gate cannot represent at rest.
        success *= initializationRejectedAtomically(makeData(),
                                                    1.2,
                                                    {0.0, 77.0, 77.0, 0.5},
                                                    "active high-value gate");

        // A non-finite field-voltage seed is rejected before any signal is
        // published.
        Fixture<ScalarT> nonfinite(makeData());
        nonfinite.attachAllInputs(77.0);
        success *= nonfinite.prepare(std::numeric_limits<RealT>::quiet_NaN());
        success *= (nonfinite.esdc1a.initialize() != 0);
        success *= scalarMatches(nonfinite.input(E::VREF), 77.0, "rejected vref preservation");

        // An invalid configuration is rejected before any state is written.
        auto invalid_data                   = makeData();
        invalid_data.parameters[Params::Ka] = 0.0;
        Fixture<ScalarT> invalid_fixture(invalid_data);
        invalid_fixture.attachAllInputs();
        success *= (invalid_fixture.esdc1a.allocate() == 0);
        poisonState(invalid_fixture, 1.2);
        const auto invalid_y  = copyVector(invalid_fixture.esdc1a.y());
        const auto invalid_yp = copyVector(invalid_fixture.esdc1a.yp());
        if (invalid_fixture.esdc1a.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: invalid configuration\n";
          success = false;
        }
        success *= vectorUnchanged(invalid_fixture.esdc1a.y(), invalid_y, "state");
        success *= vectorUnchanged(invalid_fixture.esdc1a.yp(), invalid_yp, "derivative");

        // A depressed speed input rescales the seed without rejection.
        auto admissible_speed                       = makeData();
        admissible_speed.parameters[Params::Spdmlt] = true;
        Fixture<ScalarT> speed_fixture(admissible_speed);
        speed_fixture.attachAllInputs();
        speed_fixture.input(E::OMEGA)  = -0.5;
        success                       *= speed_fixture.initialize(1.2);
        success                       *= (speed_fixture.evaluate() == 0);
        success                       *= allResidualsZero(speed_fixture.esdc1a);
        success                       *= scalarMatches(speed_fixture.esdc1a.y().getData()[I::EFDP],
                                 2.4,
                                 "rescaled EFDP");

        // The gate stays representable arbitrarily close to its operating
        // point: the ramp inverse seeds a 0.003 margin exactly.
        Fixture<ScalarT> near_gate(makeData());
        near_gate.attachAllInputs();
        success *= near_gate.initialize(1.2);
        success *= (near_gate.evaluate() == 0);
        success *= allResidualsZero(near_gate.esdc1a);

        // Summing-junction routing removes the gate constraint entirely.
        auto junction_data                    = makeData();
        junction_data.parameters[Params::UEL] = static_cast<IdxT>(2);
        Fixture<ScalarT> junction(junction_data);
        junction.attachAllInputs();
        junction.input(E::VUEL)  = 0.7;
        success                 *= junction.initialize(1.2);
        success                 *= (junction.evaluate() == 0);
        success                 *= allResidualsZero(junction.esdc1a);

        return success.report(__func__);
      }

      /// A fixed numerical answer key for all 11 ESDC1A residual rows. The
      /// expected values are literals, not a second implementation of ESDC1A.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.initialize(1.2);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.esdc1a);
        success *= (fixture.evaluate() == 0);

        // Values are pinned after an independent one-time evaluation of the
        // documented equations at setAnswerKeyState()/setAnswerKeyInputs().
        const std::array<Row, I::MAXIMUM> expected{{
            {I::EFDP, 0.04000000000000004},
            {I::VC, 0.19442890089805262},
            {I::VR, 0.13666666666666663},
            {I::VF, -0.022222222222222213},
            {I::XLL, 0.024999999999999994},
            {I::EV, -0.27999999999999986},
            {I::VLL, -0.007500000000000031},
            {I::VHV, 0.31535073999664665},
            {I::SE, -0.07541019662496842},
            {I::VFE, 0.1600000000000001},
            {I::EFD, 0.9100000000000001},
        }};

        success *= (static_cast<size_t>(fixture.esdc1a.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.esdc1a, expected);

        return success.report(__func__);
      }

      /// The transducer, summing junction, lead-lag, stabilizing feedback,
      /// and regulator anti-windup behavior at driven states with literal
      /// expectations.
      TestOutcome voltageRegulation()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeData());
        fixture.attachAllInputs();
        success *= fixture.initialize(1.2);

        // Transducer: the sensed voltage relaxes toward the bus magnitude.
        setState(fixture.esdc1a, {{I::VC, 1.1}});
        setDerivative(fixture.esdc1a, {{I::VC, 0.2}});
        success *= (fixture.evaluate() == 0);
        success *= residualsMatch(fixture.esdc1a, {{I::VC, -5.2}}, "voltage transducer");

        // The field-voltage state and the stabilizing feedback share the
        // (VR - VFE) drive.
        setState(fixture.esdc1a, {{I::VR, 0.6}, {I::VFE, 0.2}, {I::VF, 0.1}});
        setDerivative(fixture.esdc1a, {{I::EFDP, 0.1}, {I::VF, 0.05}});
        success *= (fixture.evaluate() == 0);
        success *= residualsMatch(fixture.esdc1a,
                                  {{I::EFDP, 0.7}, {I::VF, -0.13571428571428573}},
                                  "field-voltage and feedback drive");

        // Summing junction: UEL < 2 excludes the UEL input from the error.
        Fixture<ScalarT> summing(makeData());
        summing.attachAllInputs();
        success                *= summing.initialize(1.2);
        summing.input(E::VREF)  = 1.1;
        summing.input(E::VS)    = 0.05;
        summing.input(E::VUEL)  = 0.2;
        setState(summing.esdc1a, {{I::VC, 0.9}, {I::VF, 0.02}, {I::EV, 0.1}});
        success *= (summing.evaluate() == 0);
        success *= residualsMatch(summing.esdc1a, {{I::EV, 0.13}}, "summing junction");

        // UEL >= 2 routes the UEL input through the summing junction and
        // turns the high-value gate into a lead-lag passthrough.
        auto junction_data                    = makeData();
        junction_data.parameters[Params::UEL] = static_cast<IdxT>(2);
        Fixture<ScalarT> junction(junction_data);
        junction.attachAllInputs();
        success                 *= junction.initialize(1.2);
        junction.input(E::VREF)  = 1.1;
        junction.input(E::VS)    = 0.05;
        junction.input(E::VUEL)  = 0.2;
        setState(junction.esdc1a,
                 {{I::VC, 0.9}, {I::VF, 0.02}, {I::EV, 0.1}, {I::VLL, 0.5}, {I::VHV, 0.2}});
        success *= (junction.evaluate() == 0);
        success *= residualsMatch(junction.esdc1a,
                                  {{I::EV, 0.33}, {I::VHV, 0.3}},
                                  "summing-junction UEL routing");

        // An active lead-lag pair advances the error and relaxes its state.
        auto lead_lag_data                   = makeData();
        lead_lag_data.parameters[Params::Tc] = 0.2;
        Fixture<ScalarT> lead_lag(lead_lag_data);
        lead_lag.attachAllInputs();
        success *= lead_lag.initialize(1.2);
        setState(lead_lag.esdc1a, {{I::XLL, 0.4}, {I::EV, 0.7}, {I::VLL, 0.5}});
        setDerivative(lead_lag.esdc1a, {{I::XLL, 0.0}});
        success *= (lead_lag.evaluate() == 0);
        success *= residualsMatch(lead_lag.esdc1a,
                                  {{I::XLL, 0.6}, {I::VLL, 0.02}},
                                  "lead-lag");

        // The regulator anti-windup blocks outward rates at both limits and
        // admits restoring rates.
        struct AntiWindupCase
        {
          const char* label;
          RealT       vr;
          RealT       vhv;
          RealT       expected;
        };

        const std::array<AntiWindupCase, 4> antiwindup_cases{{
            {"Vrmax blocks an outward regulator rate", 1.5, 0.05, 0.0},
            {"Vrmax admits a restoring regulator rate", 1.5, 0.025, -5.0},
            {"Vrmin blocks an outward regulator rate", -1.5, -0.05, 0.0},
            {"Vrmin admits a restoring regulator rate", -1.5, -0.025, 5.0},
        }};

        for (const auto& test_case : antiwindup_cases)
        {
          setState(fixture.esdc1a, {{I::VR, test_case.vr}, {I::VHV, test_case.vhv}});
          setDerivative(fixture.esdc1a, {{I::VR, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.esdc1a,
                                    {{I::VR, test_case.expected}},
                                    test_case.label);
        }

        return success.report(__func__);
      }

      /// High-value gate selection, quadratic saturation, the exciter
      /// feedback lower limit, and the speed multiplier at driven states
      /// with literal expectations.
      TestOutcome excitationLimits()
      {
        TestStatus success = true;

        // The gate passes the larger of VLL and VUEL when UEL < 2; the
        // smooth maximum keeps two-sided sensitivity at the tie.
        struct GateCase
        {
          const char* label;
          RealT       vuel;
          RealT       expected;
        };

        const std::array<GateCase, 3> gate_cases{{
            {"gate selects the lead-lag branch", -0.5, 0.3},
            {"gate selects the UEL branch", 0.8, 0.6000000000000001},
            {"gate tie point", 0.5, 0.30288811325233306},
        }};

        Fixture<ScalarT> gate(makeData());
        gate.attachAllInputs();
        success *= gate.initialize(1.2);
        for (const auto& test_case : gate_cases)
        {
          gate.input(E::VUEL) = test_case.vuel;
          setState(gate.esdc1a, {{I::VLL, 0.5}, {I::VHV, 0.2}});
          success *= (gate.evaluate() == 0);
          success *= residualsMatch(gate.esdc1a, {{I::VHV, test_case.expected}}, test_case.label);
        }

        // Quadratic saturation above and below the fitted knee, then with
        // the fit disabled at the same field voltage.
        Fixture<ScalarT> saturation(makeResidualData());
        saturation.attachAllInputs();
        success *= saturation.initialize(1.2);
        setState(saturation.esdc1a, {{I::EFDP, 2.0}, {I::SE, 0.05}});
        success *= (saturation.evaluate() == 0);
        success *= residualsMatch(saturation.esdc1a,
                                  {{I::SE, -0.035410196624968436}},
                                  "saturation above the knee");
        setState(saturation.esdc1a, {{I::EFDP, 1.0}});
        success *= (saturation.evaluate() == 0);
        success *= residualsMatch(saturation.esdc1a,
                                  {{I::SE, -0.05}},
                                  "saturation below the knee");

        auto disabled_data                    = makeResidualData();
        disabled_data.parameters[Params::Se1] = 0.0;
        disabled_data.parameters[Params::Se2] = 0.0;
        Fixture<ScalarT> disabled(disabled_data);
        disabled.attachAllInputs();
        success *= disabled.initialize(1.2);
        setState(disabled.esdc1a, {{I::EFDP, 2.0}, {I::SE, 0.05}});
        success *= (disabled.evaluate() == 0);
        success *= residualsMatch(disabled.esdc1a, {{I::SE, -0.05}}, "saturation disabled");

        // The exciter feedback lower limit clamps a negative feedback drive
        // to zero only while enabled.
        auto feedback_data                    = makeData();
        feedback_data.parameters[Params::Ke]  = -0.2;
        feedback_data.parameters[Params::Se1] = 0.0;
        feedback_data.parameters[Params::Se2] = 0.0;

        for (const auto& [limited, expected] : std::array<std::pair<bool, RealT>, 2>{{
                 {true, 0.0},
                 {false, -0.2},
             }})
        {
          auto data                       = feedback_data;
          data.parameters[Params::exclim] = limited;
          Fixture<ScalarT> feedback(data);
          feedback.attachAllInputs();
          feedback.input(E::VUEL)  = -0.5;
          success                 *= feedback.initialize(1.2);
          setState(feedback.esdc1a, {{I::EFDP, 1.0}, {I::SE, 0.0}, {I::VFE, 0.0}});
          success *= (feedback.evaluate() == 0);
          success *= residualsMatch(feedback.esdc1a,
                                    {{I::VFE, expected}},
                                    limited ? "feedback lower limit engaged"
                                            : "feedback lower limit disabled");
        }

        // The speed multiplier scales the published field voltage only when
        // enabled.
        for (const auto& [enabled, expected] : std::array<std::pair<bool, RealT>, 2>{{
                 {false, 0.0},
                 {true, 0.06000000000000005},
             }})
        {
          auto data                       = makeData();
          data.parameters[Params::Spdmlt] = enabled;
          Fixture<ScalarT> speed(data);
          speed.attachAllInputs();
          success               *= speed.initialize(1.2);
          speed.input(E::OMEGA)  = 0.05;
          setState(speed.esdc1a, {{I::EFDP, 1.2}, {I::EFD, 1.2}});
          success *= (speed.evaluate() == 0);
          success *= residualsMatch(speed.esdc1a,
                                    {{I::EFD, expected}},
                                    enabled ? "speed multiplier enabled"
                                            : "speed multiplier disabled");
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// A single rich state and all four external inputs drive both
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
            std::cout << "ESDC1A Jacobian row " << row
                      << " mismatch between dependency tracking and Enzyme\n";
            success = false;
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Exciter::Esdc1aParameters;
      using Vars   = PhasorDynamics::Exciter::Esdc1aInternalVariables;
      using Ext    = PhasorDynamics::Exciter::Esdc1aExternalVariables;
      using Mon    = PhasorDynamics::Exciter::Esdc1aMonitorableVariables;
      using Data   = PhasorDynamics::Exciter::Esdc1aData<RealT, IdxT>;
      using I      = PhasorDynamics::Exciter::Esdc1aIdx;
      using E      = PhasorDynamics::Exciter::Esdc1aExt;

      /// A vector row paired with a value: either an input to write or an
      /// expected result. Rows are `Esdc1aIdx`/`Esdc1aExt` constants, so a
      /// failure report locates itself without any name string to maintain.
      using Row     = std::pair<size_t, RealT>;
      using Rows    = std::initializer_list<Row>;
      using Esdc1aT = PhasorDynamics::Exciter::Esdc1a<ScalarT, IdxT>;

      /// Owns the terminal bus, ESDC1A, the assigned field-voltage node, and
      /// the attached input nodes. Signal storage is declared before the
      /// model so every referenced node outlives ESDC1A. Copying would
      /// invalidate the model and signal-node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, E::MAXIMUM>                                   input_values_{};
        std::array<IdxT, E::MAXIMUM>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, E::MAXIMUM> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> efd_node_;

      public:
        explicit Fixture(const Data& data, RealT vr = 0.8, RealT vi = 0.6)
          : bus(static_cast<T>(vr), static_cast<T>(vi)),
            esdc1a(&bus, data)
        {
          esdc1a.getSignals().template assignSignalNode<Vars::EFD>(&efd_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Attach fixture-owned storage to every external input.
        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = esdc1a.size() + bus.size();

          for (size_t port = 0; port < E::MAXIMUM; ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = esdc1a.getSignals();
          signals.template attachSignalNode<Ext::OMEGA>(&input_nodes_[E::OMEGA]);
          signals.template attachSignalNode<Ext::VREF>(&input_nodes_[E::VREF]);
          signals.template attachSignalNode<Ext::VS>(&input_nodes_[E::VS]);
          signals.template attachSignalNode<Ext::VUEL>(&input_nodes_[E::VUEL]);
        }

        /// Seed the assigned field-voltage node.
        void seedEfd(RealT efd)
        {
          efd_node_.init(static_cast<T>(efd));
        }

        /// Everything ESDC1A initialization requires: allocation,
        /// verification, an initialized terminal bus, and a seeded
        /// field-voltage node.
        bool prepare(RealT efd)
        {
          const bool success = (bus.allocate() == 0) && (esdc1a.allocate() == 0)
                               && (esdc1a.verify() == 0) && (bus.initialize() == 0);
          if (!success)
          {
            std::cout << "ESDC1A fixture preparation failed\n";
            return false;
          }

          seedEfd(efd);
          return true;
        }

        /// prepare() plus successful ESDC1A initialization.
        bool initialize(RealT efd)
        {
          if (!prepare(efd))
          {
            return false;
          }
          if (esdc1a.initialize() != 0)
          {
            std::cout << "ESDC1A initialization failed\n";
            return false;
          }
          return true;
        }

        int evaluate()
        {
          return esdc1a.evaluateResidual();
        }

        T efd() const
        {
          return efd_node_.read();
        }

        T& input(size_t port)
        {
          return input_values_[port];
        }

        IdxT inputIndex(size_t port) const
        {
          return input_indices_[port];
        }

        PhasorDynamics::Bus<T, IdxT>             bus;
        PhasorDynamics::Exciter::Esdc1a<T, IdxT> esdc1a;
      };

      static constexpr RealT kStateVr = 0.9;
      static constexpr RealT kStateVi = 0.4;

      Data makeMinimalData() const
      {
        Data data;
        data.device_class          = "Esdc1a";
        data.disambiguation_string = "esdc1a_test";
        data.monitored_variables.insert(Mon::efd);
        data.monitored_variables.insert(Mon::vc);
        data.monitored_variables.insert(Mon::vr);
        data.monitored_variables.insert(Mon::vf);
        data.monitored_variables.insert(Mon::se);
        data.monitored_variables.insert(Mon::vfe);
        return data;
      }

      Data makeExplicitDefaultData() const
      {
        auto data = makeMinimalData();

        // These are the documented defaults.
        data.parameters[Params::Tr]     = 0.0;
        data.parameters[Params::Ka]     = 40.0;
        data.parameters[Params::Ta]     = 0.1;
        data.parameters[Params::Tb]     = 0.0;
        data.parameters[Params::Tc]     = 0.0;
        data.parameters[Params::Vrmax]  = 1.0;
        data.parameters[Params::Vrmin]  = -1.0;
        data.parameters[Params::Ke]     = 0.1;
        data.parameters[Params::Te]     = 0.5;
        data.parameters[Params::Kf]     = 0.05;
        data.parameters[Params::Tf1]    = 0.7;
        data.parameters[Params::Spdmlt] = false;
        data.parameters[Params::E1]     = 2.8;
        data.parameters[Params::Se1]    = 0.08;
        data.parameters[Params::E2]     = 3.7;
        data.parameters[Params::Se2]    = 0.33;
        data.parameters[Params::UEL]    = static_cast<IdxT>(0);
        data.parameters[Params::exclim] = true;
        return data;
      }

      Data makeData() const
      {
        auto data = makeMinimalData();

        // The documented typical values with the floored time constants
        // raised above the floor, so routine fixtures log no warnings.
        data.parameters[Params::Tr]     = 0.02;
        data.parameters[Params::Ka]     = 40.0;
        data.parameters[Params::Ta]     = 0.1;
        data.parameters[Params::Tb]     = 0.5;
        data.parameters[Params::Tc]     = 0.0;
        data.parameters[Params::Vrmax]  = 1.0;
        data.parameters[Params::Vrmin]  = -1.0;
        data.parameters[Params::Ke]     = 0.1;
        data.parameters[Params::Te]     = 0.5;
        data.parameters[Params::Kf]     = 0.05;
        data.parameters[Params::Tf1]    = 0.7;
        data.parameters[Params::Spdmlt] = false;
        data.parameters[Params::E1]     = 2.8;
        data.parameters[Params::Se1]    = 0.08;
        data.parameters[Params::E2]     = 3.7;
        data.parameters[Params::Se2]    = 0.33;
        data.parameters[Params::UEL]    = static_cast<IdxT>(0);
        data.parameters[Params::exclim] = true;
        return data;
      }

      Data makeResidualData() const
      {
        auto data = makeData();

        // Dynamic-response parameters: every gain, lag, and saturation
        // coefficient is nontrivial, the lead-lag is active, and the speed
        // multiplier is enabled.
        data.parameters[Params::Tr]     = 0.2;
        data.parameters[Params::Ka]     = 25.0;
        data.parameters[Params::Ta]     = 0.3;
        data.parameters[Params::Tb]     = 0.8;
        data.parameters[Params::Tc]     = 0.3;
        data.parameters[Params::Vrmax]  = 2.0;
        data.parameters[Params::Vrmin]  = -2.0;
        data.parameters[Params::Ke]     = 0.2;
        data.parameters[Params::Te]     = 0.6;
        data.parameters[Params::Kf]     = 0.08;
        data.parameters[Params::Tf1]    = 0.9;
        data.parameters[Params::Spdmlt] = true;
        data.parameters[Params::E1]     = 2.4;
        data.parameters[Params::Se1]    = 0.1;
        data.parameters[Params::E2]     = 3.2;
        data.parameters[Params::Se2]    = 0.5;
        return data;
      }

      /// The external inputs the residual answer key is evaluated against.
      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(E::OMEGA) = 0.03;
        fixture.input(E::VREF)  = 1.05;
        fixture.input(E::VS)    = 0.04;
        fixture.input(E::VUEL)  = 0.334;
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass,
      /// and VLL sits close enough to VUEL that the smooth gate keeps
      /// two-sided sensitivity.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Exciter::Esdc1a<T, IdxT>& esdc1a) const
      {
        setState(esdc1a,
                 {{I::EFDP, 2.00}, {I::VC, 0.95}, {I::VR, 0.45}, {I::VF, 0.06}, {I::XLL, 0.30}, {I::EV, 0.36}, {I::VLL, 0.33}, {I::VHV, 0.02}, {I::SE, 0.09}, {I::VFE, 0.42}, {I::EFD, 1.15}});
        setDerivative(esdc1a,
                      {{I::EFDP, 0.01},
                       {I::VC, -0.02},
                       {I::VR, 0.03},
                       {I::VF, -0.04},
                       {I::XLL, 0.05}});
      }

      /// Omitting every parameter must give exactly the model built from the
      /// defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), kStateVr, kStateVi);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), kStateVr, kStateVi);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        bool success = implicit_defaults.initialize(1.2)
                       && explicit_defaults.initialize(1.2);
        if (!success)
        {
          std::cout << "ESDC1A documented-default comparison failed to initialize\n";
          return false;
        }

        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.esdc1a.y(),
                                   copyVector(explicit_defaults.esdc1a.y()),
                                   "documented-default state");
        success *= vectorUnchanged(implicit_defaults.esdc1a.yp(),
                                   copyVector(explicit_defaults.esdc1a.yp()),
                                   "documented-default derivative");
        success *= vectorUnchanged(implicit_defaults.esdc1a.getResidual(),
                                   copyVector(explicit_defaults.esdc1a.getResidual()),
                                   "documented-default residual");

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.esdc1a);
        setAnswerKeyState(explicit_defaults.esdc1a);
        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.esdc1a.getResidual(),
                                   copyVector(explicit_defaults.esdc1a.getResidual()),
                                   "documented-default dynamic residual");
        return success;
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        Fixture<ScalarT> fixture(data);
        return fixture.esdc1a.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeData());
        fixture.esdc1a.getSignals().template attachSignalNode<variable>(&unlinked_node);
        return fixture.esdc1a.verify() > 0;
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
      /// the aliased efd entry, so any write by a rejected initialization
      /// is visible.
      void poisonState(Fixture<ScalarT>& fixture, RealT efd) const
      {
        auto* y  = fixture.esdc1a.y().getData();
        auto* yp = fixture.esdc1a.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.esdc1a.y().getSize()); ++i)
        {
          y[i]  = 0.125 + 0.01 * static_cast<RealT>(i);
          yp[i] = -0.25 - 0.01 * static_cast<RealT>(i);
        }
        fixture.seedEfd(efd);
        fixture.esdc1a.y().setDataUpdated();
        fixture.esdc1a.yp().setDataUpdated();
      }

      bool initializationRejectedAtomically(const Data&                          data,
                                            RealT                                efd_seed,
                                            const std::array<RealT, E::MAXIMUM>& inputs,
                                            const char*                          label) const
      {
        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          fixture.input(port) = inputs[port];
        }
        if (!fixture.prepare(efd_seed))
        {
          return false;
        }

        poisonState(fixture, efd_seed);
        const auto y_before  = copyVector(fixture.esdc1a.y());
        const auto yp_before = copyVector(fixture.esdc1a.yp());

        bool success = true;
        if (fixture.esdc1a.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: " << label << "\n";
          success = false;
        }

        success *= scalarMatches(fixture.efd(), efd_seed, "rejected efd preservation");
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          success &= rowMatches(static_cast<RealT>(fixture.input(port)),
                                inputs[port],
                                "external input",
                                port,
                                "changed");
        }
        success *= vectorUnchanged(fixture.esdc1a.y(), y_before, "state");
        success *= vectorUnchanged(fixture.esdc1a.yp(), yp_before, "derivative");
        return success;
      }

      /// Write state rows and publish the update, folding in the
      /// setDataUpdated() that a hand-written write block has to remember.
      template <typename T>
      void setState(PhasorDynamics::Exciter::Esdc1a<T, IdxT>& esdc1a, Rows rows) const
      {
        auto* y = esdc1a.y().getData();
        for (const auto& [row, value] : rows)
        {
          y[row] = static_cast<T>(value);
        }
        esdc1a.y().setDataUpdated();
      }

      /// setState() for the derivative vector.
      template <typename T>
      void setDerivative(PhasorDynamics::Exciter::Esdc1a<T, IdxT>& esdc1a, Rows rows) const
      {
        auto* yp = esdc1a.yp().getData();
        for (const auto& [row, value] : rows)
        {
          yp[row] = static_cast<T>(value);
        }
        esdc1a.yp().setDataUpdated();
      }

      /// Compare one vector row against its expected value. Every row check
      /// in this suite reports through here, so failures share one format.
      /// Rows are named by position, which is the `Esdc1aIdx` constant the
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
        std::cout << "ESDC1A " << what << " row " << row << ' ' << context
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

      bool residualsMatch(const Esdc1aT& esdc1a, Rows rows, const char* context = "") const
      {
        return rowsMatch(esdc1a.getResidual(), rows.begin(), rows.size(), "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const Esdc1aT&               esdc1a,
                          const std::array<Row, size>& rows,
                          const char*                  context = "") const
      {
        return rowsMatch(esdc1a.getResidual(), rows.data(), size, "residual", context);
      }

      /// The model sits at a steady state: every residual and every
      /// derivative is zero.
      bool allResidualsZero(const Esdc1aT& esdc1a) const
      {
        bool        success = true;
        const auto* f       = esdc1a.getResidual().getData();
        const auto* yp      = esdc1a.yp().getData();
        for (size_t row = 0; row < static_cast<size_t>(esdc1a.getResidual().getSize()); ++row)
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
        auto* y     = fixture.esdc1a.y().getData();
        auto* yp    = fixture.esdc1a.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        const auto model_size = static_cast<size_t>(fixture.esdc1a.size());
        for (size_t i = 0; i < model_size; ++i)
        {
          y[i].setVariableNumber(i);
          yp[i].setVariableNumber(i);
        }
        for (size_t i = 0; i < static_cast<size_t>(fixture.bus.size()); ++i)
        {
          bus_y[i].setVariableNumber(model_size + i);
        }
        for (size_t port = 0; port < E::MAXIMUM; ++port)
        {
          fixture.input(port).setVariableNumber(fixture.inputIndex(port));
        }

        fixture.esdc1a.y().setDataUpdated();
        fixture.esdc1a.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          const Data& data,
          TestStatus& success) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.initialize(1.2);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.esdc1a);
        numberVariables(fixture);
        success *= (fixture.evaluate() == 0);

        const auto                         model_size = static_cast<size_t>(fixture.esdc1a.size());
        std::vector<DepVar::DependencyMap> rows(model_size);
        const auto*                        f = fixture.esdc1a.getResidual().getData();
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
        Fixture<ScalarT> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.initialize(1.2);

        for (IdxT i = 0; i < fixture.bus.size(); ++i)
        {
          fixture.bus.setVariableIndex(i, fixture.esdc1a.size() + i);
        }

        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.esdc1a);
        fixture.esdc1a.updateTime(0.0, 1.0);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.esdc1a.evaluateJacobian() == 0);
        success *= (fixture.esdc1a.constructCsr() == 0);
        return MapFromCsr(fixture.esdc1a.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
