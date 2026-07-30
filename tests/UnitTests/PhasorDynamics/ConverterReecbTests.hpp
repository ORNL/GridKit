#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/ReecbData.hpp>
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
    class ConverterReecbTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterReecbTests()  = default;
      ~ConverterReecbTests() = default;

      // Fixed answer keys should agree to ordinary double-precision roundoff.
      static constexpr RealT kBehaviorTol = 1.0e-12;

      // REECB initialization solves smooth-limiter inputs near their asymptotes.
      // The resulting steady residuals are O(1e-10).
      static constexpr RealT kSteadyStateTol = 1.0e-9;

      // Enzyme and dependency tracking traverse the same smooth expressions
      // differently; their double-precision derivatives agree to O(1e-10).
      static constexpr RealT kJacobianTol = 1.0e-9;

      /// Construction and every verify() error class, including parameter
      /// types, parameter relationships, bus ownership, and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(index(Vars::MAXIMUM)));
        success *= (empty.getMonitor() == nullptr);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> configured(&bus, makeData());
        success *= (configured.size() == static_cast<IdxT>(index(Vars::MAXIMUM)));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() == 0);

        noteExpectedLogs("Testing REECB defaults and invalid configurations. "
                         "Logged errors and time-constant warnings are expected.");

        auto minimal_data                    = makeMinimalData();
        minimal_data.parameters[Params::mva] = 100.0;
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> minimal(&bus, minimal_data);
        success *= (minimal.verify() == 0);
        success *= defaultsMatchDocumentedValues();

        success *= (empty.verify() > 0);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> missing_mva(&bus, makeMinimalData());
        success *= (missing_mva.verify() > 0);

        success *= invalidParameterCase(bus, Params::mva, 0.0);
        success *= invalidParameterCase(bus, Params::Trv, -0.1);
        success *= invalidParameterCase(bus, Params::Tp, -0.1);
        success *= invalidParameterCase(bus, Params::Tiq, -0.1);
        success *= invalidParameterCase(bus, Params::Tpord, -0.1);
        success *= invalidParameterCase(bus, Params::Vdip, 1.2);
        success *= invalidParameterCase(bus, Params::dbd1, 0.1);
        success *= invalidParameterCase(bus, Params::dbd2, -0.1);
        success *= invalidParameterCase(bus, Params::Iql1, 2.0);
        success *= invalidParameterCase(bus, Params::Qmin, 2.0);
        success *= invalidParameterCase(bus, Params::Vmin, 2.0);
        success *= invalidParameterCase(bus, Params::dPmin, 0.0);
        success *= invalidParameterCase(bus, Params::dPmax, 0.0);
        success *= invalidParameterCase(bus, Params::Pmin, 2.0);
        success *= invalidParameterCase(bus, Params::Imax, -0.1);

        for (const Params flag : {Params::PfFlag, Params::VFlag, Params::QFlag, Params::Pqflag})
        {
          auto bad_integer             = makeData();
          bad_integer.parameters[flag] = static_cast<IdxT>(2);
          PhasorDynamics::Converter::Reecb<ScalarT, IdxT> bad_integer_model(&bus, bad_integer);
          success *= (bad_integer_model.verify() > 0);

          // Real-valued 0/1 is intentionally rejected: switches are JSON
          // booleans or integer 0/1, matching REGCA's parameter contract.
          auto bad_real             = makeData();
          bad_real.parameters[flag] = static_cast<RealT>(1.0);
          PhasorDynamics::Converter::Reecb<ScalarT, IdxT> bad_real_model(&bus, bad_real);
          success *= (bad_real_model.verify() > 0);
        }

        auto integer_switches                       = makeData();
        integer_switches.parameters[Params::PfFlag] = static_cast<IdxT>(0);
        integer_switches.parameters[Params::VFlag]  = static_cast<IdxT>(1);
        integer_switches.parameters[Params::QFlag]  = static_cast<IdxT>(0);
        integer_switches.parameters[Params::Pqflag] = static_cast<IdxT>(1);
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> integer_switch_model(
            &bus,
            integer_switches);
        success *= (integer_switch_model.verify() == 0);

        auto bad_numeric_type                    = makeData();
        bad_numeric_type.parameters[Params::mva] = true;
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> bad_numeric_model(
            &bus,
            bad_numeric_type);
        success *= (bad_numeric_model.verify() > 0);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> busless(nullptr, makeData());
        success *= (busless.verify() > 0);

        success *= unlinkedSignalRejected<Ext::PE>(bus);
        success *= unlinkedSignalRejected<Ext::QGEN>(bus);
        success *= unlinkedSignalRejected<Ext::QEXT>(bus);
        success *= unlinkedSignalRejected<Ext::PFAREF>(bus);
        success *= unlinkedSignalRejected<Ext::PREF>(bus);

        // All four zero time constants use the documented numerical floor and
        // still admit a consistent steady-state initialization.
        auto zero_time                      = makeData();
        zero_time.parameters[Params::Trv]   = 0.0;
        zero_time.parameters[Params::Tp]    = 0.0;
        zero_time.parameters[Params::Tiq]   = 0.0;
        zero_time.parameters[Params::Tpord] = 0.0;

        Fixture<ScalarT> fixture(zero_time);
        success *= fixture.initialize(0.2, 0.6);
        success *= (fixture.evaluate() == 0);
        success *= allResidualsZero(fixture.reecb);

        return success.report(__func__);
      }

      /// A nonidentity power-base initialization with every port attached.
      /// Assigned command nodes are seeded after allocate() and must remain
      /// unchanged while REECB initializes and publishes its feedback signals.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        auto data                    = makeData();
        data.parameters[Params::mva] = 50.0;

        Fixture<ScalarT> fixture(data, 0.8, 0.6);
        fixture.attachAllInputs(99.0);
        success *= fixture.initialize(0.05, 0.25);
        success *= (fixture.reecb.tagDifferentiable() == 0);
        success *= (fixture.evaluate() == 0);

        const auto* y  = fixture.reecb.y().getData();
        success       *= scalarMatches(y[index(Vars::VT)], 1.0, "VT");
        success       *= scalarMatches(y[index(Vars::VMEAS)], 1.0, "VMEAS");
        success       *= scalarMatches(y[index(Vars::PMEAS)], 0.5, "PMEAS on component base");
        success       *= scalarMatches(y[index(Vars::QREF)], 0.1, "QREF on component base");
        success       *= scalarMatches(y[index(Vars::PORD)], 0.5, "PORD on component base");
        success       *= scalarMatches(fixture.iqcmd(), 0.05, "seeded iqcmd");
        success       *= scalarMatches(fixture.ipcmd(), 0.25, "seeded ipcmd");

        success *= scalarMatches(fixture.input(Ext::PE), 0.25, "published pe");
        success *= scalarMatches(fixture.input(Ext::QGEN), 0.05, "published qgen");
        success *= scalarMatches(fixture.input(Ext::QEXT), 0.05, "published qext");
        success *= scalarMatches(fixture.input(Ext::PFAREF), 0.0, "inactive pfaref fallback");
        success *= scalarMatches(fixture.input(Ext::PREF), 0.25, "published pref");

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.reecb.getMonitor());
        std::stringstream monitor_output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, monitor_output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string monitor_header;
        std::string monitor_values;
        std::getline(monitor_output, monitor_header);
        std::getline(monitor_output, monitor_values);
        success              *= (monitor_header == "t,Reecb_reecb_test_iqcmd,Reecb_reecb_test_ipcmd,"
                                                   "Reecb_reecb_test_vmeas,Reecb_reecb_test_pmeas");
        const auto monitored  = Tokenizer<RealT>(monitor_values, ',')();
        if (monitored.size() == 5)
        {
          success *= scalarMatches(monitored[1], 0.05, "monitored iqcmd");
          success *= scalarMatches(monitored[2], 0.25, "monitored ipcmd");
          success *= scalarMatches(monitored[3], 1.0, "monitored vmeas");
          success *= scalarMatches(monitored[4], 0.5, "monitored pmeas");
        }
        else
        {
          std::cout << "REECB monitor emitted " << monitored.size()
                    << " values instead of 5\n";
          success = false;
        }

        for (size_t i = 0; i < static_cast<size_t>(fixture.reecb.size()); ++i)
        {
          const bool expected = i <= index(Vars::PORD);
          if (fixture.reecb.tag()[i] != expected)
          {
            std::cout << "REECB differentiability tag " << i << " mismatch\n";
            success = false;
          }
        }
        success *= allResidualsZero(fixture.reecb);

        // Every flag combination must preserve the seeded commands and produce
        // a zero-derivative, zero-residual state. The nonzero PI gains ensure
        // that active and parked anti-windup paths are both initialized.
        for (const bool pf_flag : {false, true})
        {
          for (const bool v_flag : {false, true})
          {
            for (const bool q_flag : {false, true})
            {
              for (const bool p_priority : {false, true})
              {
                auto scenario_data                       = makeDynamicData();
                scenario_data.parameters[Params::PfFlag] = pf_flag;
                scenario_data.parameters[Params::VFlag]  = v_flag;
                scenario_data.parameters[Params::QFlag]  = q_flag;
                scenario_data.parameters[Params::Pqflag] = p_priority;
                scenario_data.parameters[Params::kqv]    = 0.0;

                Fixture<ScalarT> scenario(scenario_data);
                scenario.attachAllInputs(99.0);
                if (!scenario.initialize(0.05, 0.2))
                {
                  std::cout << "REECB initialization scenario failed: PfFlag=" << pf_flag
                            << ", VFlag=" << v_flag
                            << ", QFlag=" << q_flag
                            << ", Pqflag=" << p_priority << '\n';
                  success = false;
                  continue;
                }

                success *= (scenario.evaluate() == 0);
                success *= allResidualsZero(scenario.reecb);
                success *= scalarMatches(scenario.iqcmd(), 0.05, "scenario iqcmd preservation");
                success *= scalarMatches(scenario.ipcmd(), 0.2, "scenario ipcmd preservation");
                success *= scalarMatches(scenario.input(Ext::PE), 0.2, "scenario pe publication");
                success *= scalarMatches(scenario.input(Ext::QGEN), 0.05, "scenario qgen publication");
              }
            }
          }
        }

        // An in-band point and both voltage-band exits exercise the direct
        // Iq-injection compensation during initialization.
        for (const RealT terminal_voltage :
             {static_cast<RealT>(0.6), static_cast<RealT>(1.1), static_cast<RealT>(1.3)})
        {
          auto voltage_data                      = makeDynamicData();
          voltage_data.parameters[Params::QFlag] = true;
          voltage_data.parameters[Params::VFlag] = true;
          voltage_data.parameters[Params::Vref0] = 1.0;
          voltage_data.parameters[Params::Vmin]  = 0.5;
          voltage_data.parameters[Params::Vmax]  = 1.4;

          Fixture<ScalarT> voltage_scenario(voltage_data, terminal_voltage);
          voltage_scenario.attachAllInputs();
          success *= voltage_scenario.initialize(0.05, 0.2);
          success *= (voltage_scenario.evaluate() == 0);
          success *= allResidualsZero(voltage_scenario.reecb);
          success *= scalarMatches(voltage_scenario.iqcmd(), 0.05, "voltage-band iqcmd preservation");
          success *= scalarMatches(voltage_scenario.ipcmd(), 0.2, "voltage-band ipcmd preservation");
        }

        return success.report(__func__);
      }

      /// Current, power, and selected-controller initialization domains.
      /// Every rejection is atomic; zero-current and exact limiter boundaries
      /// remain admissible.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible REECB current, power, and controller "
                         "initialization points. Logged errors are expected.");

        struct RejectionCase
        {
          const char* label;
          RealT       iqcmd;
          RealT       ipcmd;
          RealT       imax;
          RealT       pmin;
          RealT       pmax;
        };

        const std::array<RejectionCase, 4> rejected{{
            {"negative active-current command", 0.0, -0.1, 1.0, 0.0, 1.0},
            {"command vector outside Imax", 0.8, 0.8, 1.0, 0.0, 1.0},
            {"active-power seed above Pmax", 0.0, 0.8, 1.0, 0.0, 0.5},
            {"active-power seed below Pmin", 0.0, 0.2, 1.0, 0.3, 1.0},
        }};

        for (const auto& test_case : rejected)
        {
          auto data                      = makeData();
          data.parameters[Params::Imax]  = test_case.imax;
          data.parameters[Params::Pmin]  = test_case.pmin;
          data.parameters[Params::Pmax]  = test_case.pmax;
          success                       *= initializationRejectedAtomically(
              data, test_case.iqcmd, test_case.ipcmd, 1.0, test_case.label);
        }

        // With reactive-current control selected, the voltage-controller
        // limiter input must reproduce the compensated reactive-current target.
        auto controller_data                       = makeData();
        controller_data.parameters[Params::QFlag]  = true;
        controller_data.parameters[Params::Imax]   = 1.0;
        controller_data.parameters[Params::Vref0]  = 1.0;
        controller_data.parameters[Params::kqv]    = 10.0;
        controller_data.parameters[Params::Iql1]   = -1.0;
        controller_data.parameters[Params::Iqh1]   = 1.0;
        success                                   *= initializationRejectedAtomically(
            controller_data,
            -0.5,
            0.0,
            0.6,
            "reactive-current target outside the selected voltage-controller limits");

        // VFlag requires limiter inputs for terminal voltage through Vmin/Vmax
        // and initial reactive power through Qmin/Qmax.
        auto voltage_data                       = makeData();
        voltage_data.parameters[Params::QFlag]  = true;
        voltage_data.parameters[Params::VFlag]  = true;
        voltage_data.parameters[Params::Vmax]   = 0.9;
        success                                *= initializationRejectedAtomically(
            voltage_data,
            0.1,
            0.2,
            1.0,
            "terminal voltage outside selected Vmin/Vmax");

        auto reactive_power_data                       = makeData();
        reactive_power_data.parameters[Params::QFlag]  = true;
        reactive_power_data.parameters[Params::VFlag]  = true;
        reactive_power_data.parameters[Params::Qmin]   = -0.1;
        reactive_power_data.parameters[Params::Qmax]   = 0.1;
        success                                       *= initializationRejectedAtomically(
            reactive_power_data,
            0.2,
            0.2,
            1.0,
            "initial reactive power outside selected Qmin/Qmax");

        // Power-factor control cannot infer an angle for nonzero Q at zero
        // active power. This late rejection proves initialization atomicity.
        auto power_factor_data                        = makeData();
        power_factor_data.parameters[Params::PfFlag]  = true;
        success                                      *= initializationRejectedAtomically(
            power_factor_data,
            0.2,
            0.0,
            1.0,
            "nonzero power-factor reactive reference at zero active power");

        // Zero current and both exact current-circle boundaries stay admissible.
        struct AdmissibleCase
        {
          RealT iqcmd;
          RealT ipcmd;
          RealT imax;
        };

        for (const auto& accepted : std::array<AdmissibleCase, 3>{{
                 {0.0, 0.0, 1.0},
                 {0.0, 1.0, 1.0},
                 {1.0, 0.0, 1.0},
             }})
        {
          auto data                     = makeData();
          data.parameters[Params::Imax] = accepted.imax;
          data.parameters[Params::Pmax] = 1.0;

          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize(accepted.iqcmd, accepted.ipcmd);
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.reecb);
        }

        return success.report(__func__);
      }

      /// A fixed numerical answer key for all 25 REECB residual rows. The
      /// expected values are literals, not a second implementation of REECB.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeDynamicData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.1, 0.2);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.reecb);
        success *= (fixture.evaluate() == 0);

        // Values are pinned after an independent one-time evaluation of the
        // documented equations at setAnswerKeyState()/setAnswerKeyInputs().
        const std::array<Row, index(Vars::MAXIMUM)> expected{{
            {index(Vars::VMEAS), 0.2400000000000002},
            {index(Vars::PMEAS), 0.1449999999999998},
            {index(Vars::XPIQ), 0.03399999970642038},
            {index(Vars::XPIV), 0.0},
            {index(Vars::QV), 0.1222222222222222},
            {index(Vars::PORD), 0.26},
            {index(Vars::VT), -0.02999999999999992},
            {index(Vars::VMEASSAFE), -0.01000000000000001},
            {index(Vars::SDIP), 0.2},
            {index(Vars::VERR), 2.821917786596795e-7},
            {index(Vars::IQV), -0.06999999999999998},
            {index(Vars::QREF), -0.2668756300679377},
            {index(Vars::EQ), 0.3499999999999999},
            {index(Vars::VPIQ), -0.1000000000000002},
            {index(Vars::EPIV), -0.04999999999999993},
            {index(Vars::FPORD), -0.02500000000000002},
            {index(Vars::RPORD), 0.05000000000000004},
            {index(Vars::IQCIRC), 0.3999999999999997},
            {index(Vars::IPCIRC), 0.8100000000000001},
            {index(Vars::IQMAX), 0.1000000000000001},
            {index(Vars::IPMAX), 0.2},
            {index(Vars::IQBASE), -0.3699999999999998},
            {index(Vars::IQRAW), -0.05000000000000002},
            {index(Vars::IQCMD), -0.1000000000000001},
            {index(Vars::IPCMD), -0.1229166666666667},
        }};

        success *= (static_cast<size_t>(fixture.reecb.getResidual().getSize()) == expected.size());
        success *= residualsMatch(fixture.reecb, expected);

        return success.report(__func__);
      }

      /// Flag selection, voltage/deadband behavior, injection limiting,
      /// reactive lag, and upper/lower/restoring anti-windup behavior.
      TestOutcome reactiveControl()
      {
        TestStatus success = true;

        struct FlagCase
        {
          const char* label;
          bool        pf;
          bool        voltage;
          bool        reactive;
          RealT       qref;
          RealT       epiv;
          RealT       iqraw;
        };

        // Toggle exactly one selector at a time so an accidental swap between
        // PfFlag, VFlag, and QFlag cannot satisfy the same answer key.
        const std::array<FlagCase, 4> cases{{
            {"all-off selectors", false, false, false, 0.4, -0.55, 0.4},
            {"PfFlag-only selectors", true, false, false, 0.11149051952976989, -0.8385094804702301, 0.4},
            {"VFlag-only selectors", false, true, false, 0.4, 0.05, 0.4},
            {"QFlag-only selectors", false, false, true, 0.4, -0.55, 0.3},
        }};

        for (const auto& test_case : cases)
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = test_case.pf;
          data.parameters[Params::VFlag]  = test_case.voltage;
          data.parameters[Params::QFlag]  = test_case.reactive;

          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          success *= fixture.initialize(0.1, 0.2);

          fixture.input(Ext::PFAREF) = 0.2;
          fixture.input(Ext::QEXT)   = 0.2; // 0.4 on the 50 MVA component base.
          setState(fixture.reecb,
                   {{index(Vars::PMEAS), 0.55},
                    {index(Vars::VMEAS), 0.95},
                    {index(Vars::VPIQ), 1.0},
                    {index(Vars::QREF), test_case.qref},
                    {index(Vars::IQBASE), 0.2},
                    {index(Vars::QV), 0.3},
                    {index(Vars::SDIP), 0.9},
                    {index(Vars::IQV), 0.1},
                    {index(Vars::EPIV), test_case.epiv},
                    {index(Vars::IQRAW), test_case.iqraw}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{index(Vars::QREF), 0.0}, {index(Vars::EPIV), 0.0}, {index(Vars::IQRAW), 0.0}},
                                    test_case.label);
        }

        Fixture<ScalarT> limit_fixture(makeDynamicData());
        limit_fixture.attachAllInputs();
        success                        *= limit_fixture.initialize(0.1, 0.2);
        limit_fixture.input(Ext::QGEN)  = 0.0;

        // A Q reference driven past each limit; EQ is the clamped reference
        // less the zeroed qgen feedback.
        for (const auto& [qref, expected_eq] : std::array<DrivenCase, 2>{{
                 {2.0, 0.8},
                 {-2.0, -0.7},
             }})
        {
          setState(limit_fixture.reecb, {{index(Vars::QREF), qref}, {index(Vars::EQ), 0.0}});
          success *= (limit_fixture.evaluate() == 0);
          success *= residualsMatch(limit_fixture.reecb,
                                    {{index(Vars::EQ), expected_eq}},
                                    "reactive-power limit");
        }

        // The same sweep through the reactive-power PI output limits.
        for (const auto& [eq, expected_vpiq] : std::array<DrivenCase, 2>{{
                 {4.0, 1.3},
                 {-4.0, 0.7},
             }})
        {
          setState(limit_fixture.reecb, {{index(Vars::EQ), eq}, {index(Vars::XPIQ), 0.0}, {index(Vars::VPIQ), 0.0}});
          success *= (limit_fixture.evaluate() == 0);
          success *= residualsMatch(limit_fixture.reecb,
                                    {{index(Vars::VPIQ), expected_vpiq}},
                                    "reactive-power PI voltage limit");
        }

        // Each row pins the smooth voltage-band/deadband/injection behavior
        // without calling CommonMath in the expected-value path.
        struct VoltageCase
        {
          RealT voltage;
          RealT vmeas;
          RealT expected_sdip_rhs;
          RealT expected_verr_rhs;
          RealT expected_iqv_rhs;
        };

        const std::array<VoltageCase, 3> voltage_cases{{
            {0.6, 0.6, 0.0, 0.37, 0.5},
            {1.0, 1.0, 1.0, -3.104066702683552e-5, -6.208133405367104e-5},
            {1.3, 1.3, 0.0, -0.28, -0.4},
        }};

        auto voltage_data                      = makeDynamicData();
        voltage_data.parameters[Params::Vref0] = 1.0;
        Fixture<ScalarT> voltage_fixture(voltage_data);
        success *= voltage_fixture.initialize(0.0, 0.2);
        for (const auto& test_case : voltage_cases)
        {
          setState(voltage_fixture.reecb,
                   {{index(Vars::VT), test_case.voltage},
                    {index(Vars::VMEAS), test_case.vmeas},
                    {index(Vars::SDIP), test_case.expected_sdip_rhs},
                    {index(Vars::VERR), test_case.expected_verr_rhs},
                    {index(Vars::IQV), test_case.expected_iqv_rhs}});
          success *= (voltage_fixture.evaluate() == 0);
          success *= residualsMatch(voltage_fixture.reecb,
                                    {{index(Vars::SDIP), 0.0}, {index(Vars::VERR), 0.0}, {index(Vars::IQV), 0.0}},
                                    "voltage band, deadband, and injection limit",
                                    kSteadyStateTol);
        }

        // The QV lag and both PI anti-windup rows are evaluated at three
        // controller directions: upper saturation, lower saturation, and a
        // restoring direction. The expected derivatives are fixed literals.
        struct AntiWindupCase
        {
          RealT xpiq;
          RealT eq;
          RealT xpiv;
          RealT epiv;
          RealT expected_xpiq;
          RealT expected_xpiv;
        };

        const std::array<AntiWindupCase, 3> antiwindup_cases{{
            {2.0, 0.5, 2.0, 0.5, 0.0, 0.0},
            {-2.0, -0.5, -2.0, -0.5, 0.0, 0.0},
            {2.0, -0.5, 2.0, -0.5, -0.2, -0.25},
        }};

        Fixture<ScalarT> controller_fixture(makeDynamicData());
        success *= controller_fixture.initialize(0.0, 0.2);
        for (const auto& test_case : antiwindup_cases)
        {
          setState(controller_fixture.reecb,
                   {{index(Vars::SDIP), 1.0},
                    {index(Vars::XPIQ), test_case.xpiq},
                    {index(Vars::EQ), test_case.eq},
                    {index(Vars::XPIV), test_case.xpiv},
                    {index(Vars::EPIV), test_case.epiv},
                    {index(Vars::IQMAX), 1.0}});
          success *= (controller_fixture.evaluate() == 0);
          success *= residualsMatch(controller_fixture.reecb,
                                    {{index(Vars::XPIQ), test_case.expected_xpiq},
                                     {index(Vars::XPIV), test_case.expected_xpiv}},
                                    "anti-windup");
        }

        setState(controller_fixture.reecb,
                 {{index(Vars::SDIP), 1.0}, {index(Vars::QREF), 0.4}, {index(Vars::VMEASSAFE), 1.0}, {index(Vars::QV), 0.2}});
        setDerivative(controller_fixture.reecb, {{index(Vars::QV), 0.1}});
        success *= (controller_fixture.evaluate() == 0);
        success *= residualsMatch(controller_fixture.reecb,
                                  {{index(Vars::QV), 0.5666666666666667}},
                                  "reactive-current lag");

        return success.report(__func__);
      }

      /// Active-power measurement/order filters, both ramp limits, power
      /// bounds, and the safe-voltage current conversion.
      TestOutcome activePowerControl()
      {
        TestStatus success = true;

        auto data                      = makeDynamicData();
        data.parameters[Params::dPmax] = 0.4;
        data.parameters[Params::dPmin] = -0.3;
        data.parameters[Params::Pmax]  = 1.0;
        data.parameters[Params::Pmin]  = 0.1;
        data.parameters[Params::Tpord] = 0.5;
        data.parameters[Params::Tp]    = 0.25;

        Fixture<ScalarT> fixture(data);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.0, 0.2);

        struct RampCase
        {
          RealT pord;
          RealT pref_system;
          RealT fpord;
          RealT rpord;
          RealT expected_fpord;
          RealT expected_rpord;
          RealT expected_pord;
        };

        const std::array<RampCase, 3> cases{{
            {0.5, 0.5, 1.0, 0.4, 0.0, 0.0, 0.4},
            {0.7, 0.1, -1.0, -0.3, 0.0, 0.0, -0.3},
            {1.2, 0.5, -0.4, -0.3, 0.0, 0.0, -0.3},
        }};

        for (const auto& test_case : cases)
        {
          fixture.input(Ext::PREF) = test_case.pref_system;
          setState(fixture.reecb,
                   {{index(Vars::PORD), test_case.pord},
                    {index(Vars::FPORD), test_case.fpord},
                    {index(Vars::RPORD), test_case.rpord},
                    {index(Vars::SDIP), 1.0}});
          setDerivative(fixture.reecb, {{index(Vars::PORD), 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{index(Vars::FPORD), test_case.expected_fpord},
                                     {index(Vars::RPORD), test_case.expected_rpord},
                                     {index(Vars::PORD), test_case.expected_pord}},
                                    "active-power order");
        }

        struct LowerPowerBoundCase
        {
          RealT       rate;
          RealT       expected_residual;
          const char* label;
        };

        const std::array<LowerPowerBoundCase, 2> lower_bound_cases{{
            {-0.3, 0.0, "Pmin blocks an outward active-power rate"},
            {0.4, 0.4, "Pmin admits a restoring active-power rate"},
        }};

        for (const auto& test_case : lower_bound_cases)
        {
          setState(fixture.reecb,
                   {{index(Vars::PORD), -0.2}, {index(Vars::RPORD), test_case.rate}, {index(Vars::SDIP), 1.0}});
          setDerivative(fixture.reecb, {{index(Vars::PORD), 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{index(Vars::PORD), test_case.expected_residual}},
                                    test_case.label);
        }

        // PE is 0.3 on system base and 0.6 on the 50 MVA component base.
        fixture.input(Ext::PE) = 0.3;
        setState(fixture.reecb,
                 {{index(Vars::PMEAS), 0.5},
                  {index(Vars::VMEAS), 0.005},
                  {index(Vars::VMEASSAFE), 0.01},
                  {index(Vars::PORD), 0.004},
                  {index(Vars::IPMAX), 1.0},
                  {index(Vars::IPCMD), 0.2}});
        setDerivative(fixture.reecb, {{index(Vars::PMEAS), 0.1}});
        success *= (fixture.evaluate() == 0);
        success *= residualsMatch(fixture.reecb,
                                  {{index(Vars::PMEAS), 0.3},
                                   {index(Vars::VMEASSAFE), 0.00109701028057513},
                                   {index(Vars::IPCMD), 0.0}},
                                  "safe-voltage active-power path");

        return success.report(__func__);
      }

      /// P- and Q-priority at the current-circle boundary. In both cases the
      /// low-priority command exactly binds the remaining-current root.
      TestOutcome currentPriority()
      {
        TestStatus success = true;

        struct PriorityCase
        {
          const char* label;
          bool        p_priority;
          RealT       iqcmd;
          RealT       ipcmd;
          RealT       iqcirc;
          RealT       ipcirc;
          RealT       iqmax;
          RealT       ipmax;
        };

        const std::array<PriorityCase, 2> cases{{
            {"Q-priority", false, 0.8, 0.6, 1.0, 0.6, 1.0, 0.6},
            {"P-priority", true, 0.6, 0.8, 0.6, 1.0, 0.6, 1.0},
        }};

        for (const auto& test_case : cases)
        {
          auto data                       = makeData();
          data.parameters[Params::Pqflag] = test_case.p_priority;
          data.parameters[Params::Imax]   = 1.0;

          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize(test_case.iqcmd, test_case.ipcmd);
          success *= (fixture.evaluate() == 0);

          success *= stateMatches(fixture.reecb,
                                  {{index(Vars::IQCIRC), test_case.iqcirc},
                                   {index(Vars::IPCIRC), test_case.ipcirc},
                                   {index(Vars::IQMAX), test_case.iqmax},
                                   {index(Vars::IPMAX), test_case.ipmax}},
                                  test_case.label);
          success *= residualsMatch(fixture.reecb,
                                    {{index(Vars::IQCMD), 0.0}, {index(Vars::IPCMD), 0.0}},
                                    test_case.label);
          success *= scalarMatches(fixture.iqcmd(), test_case.iqcmd, "priority iqcmd preservation");
          success *= scalarMatches(fixture.ipcmd(), test_case.ipcmd, "priority ipcmd preservation");
          success *= allResidualsZero(fixture.reecb);
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Representative selector modes drive both sensitivity paths. Every
      /// Enzyme CSR row must match dependency tracking, and fixed derivatives
      /// independently pin the four selector contracts.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        struct SelectorCase
        {
          const char* label;
          bool        pf;
          bool        voltage;
          bool        reactive;
          bool        p_priority;
        };

        // Start with every selector off, then activate one at a time. This keeps
        // each expected structural change attributable to exactly one flag.
        const std::array<SelectorCase, 5> cases{{
            {"all selectors off", false, false, false, false},
            {"power-factor control", true, false, false, false},
            {"voltage control", false, true, false, false},
            {"reactive-current control", false, false, true, false},
            {"active-current priority", false, false, false, true},
        }};

        for (const auto& test_case : cases)
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = test_case.pf;
          data.parameters[Params::VFlag]  = test_case.voltage;
          data.parameters[Params::QFlag]  = test_case.reactive;
          data.parameters[Params::Pqflag] = test_case.p_priority;

          const auto dependency_jacobian = dependencyTrackingJacobian(data, success);
          const auto enzyme_jacobian     = enzymeJacobian(data, success);

          success         *= (dependency_jacobian.size() == index(Vars::MAXIMUM));
          success         *= (dependency_jacobian.size() == enzyme_jacobian.size());
          const auto rows  = std::min(dependency_jacobian.size(), enzyme_jacobian.size());
          for (size_t row = 0; row < rows; ++row)
          {
            if (!isEqual(dependency_jacobian[row], enzyme_jacobian[row], kJacobianTol))
            {
              std::cout << "REECB Jacobian row " << row
                        << " mismatch between dependency tracking and Enzyme for "
                        << test_case.label << '\n';
              success = false;
            }
          }

          if (dependency_jacobian.size() != index(Vars::MAXIMUM))
          {
            continue;
          }

          // Fixed derivatives at setAnswerKeyState()/setAnswerKeyInputs(). The
          // component/system power-base ratio is two in makeDynamicData().
          RealT qext_derivative   = 2.0;
          RealT pfaref_derivative = 0.0;
          if (test_case.pf)
          {
            qext_derivative   = 0.0;
            pfaref_derivative = 0.5625630197756407;
          }

          RealT vpiq_derivative = 0.0;
          RealT qref_derivative = 1.0;
          if (test_case.voltage)
          {
            vpiq_derivative = 1.0;
            qref_derivative = 0.0;
          }

          RealT iqbase_derivative = 0.0;
          RealT qv_derivative     = 1.0;
          if (test_case.reactive)
          {
            iqbase_derivative = 1.0;
            qv_derivative     = 0.0;
          }

          RealT iqcirc_ipcmd_derivative = 0.0;
          RealT ipcirc_iqcmd_derivative = -2.0;
          if (test_case.p_priority)
          {
            iqcirc_ipcmd_derivative = -3.2;
            ipcirc_iqcmd_derivative = 0.0;
          }

          success *= dependencyMatches(dependency_jacobian,
                                       Vars::QREF,
                                       Ext::QEXT,
                                       qext_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::QREF,
                                       Ext::PFAREF,
                                       pfaref_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::EPIV,
                                       Vars::VPIQ,
                                       vpiq_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::EPIV,
                                       Vars::QREF,
                                       qref_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::IQRAW,
                                       Vars::IQBASE,
                                       iqbase_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::IQRAW,
                                       Vars::QV,
                                       qv_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::IQCIRC,
                                       Vars::IPCMD,
                                       iqcirc_ipcmd_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::IPCIRC,
                                       Vars::IQCMD,
                                       ipcirc_iqcmd_derivative,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::IQCMD,
                                       Vars::IQCMD,
                                       -2.0,
                                       test_case.label);
          success *= dependencyMatches(dependency_jacobian,
                                       Vars::IPCMD,
                                       Vars::IPCMD,
                                       -2.0,
                                       test_case.label);
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Converter::ReecbParameters;
      using Vars   = PhasorDynamics::Converter::ReecbInternalVariables;
      using Ext    = PhasorDynamics::Converter::ReecbExternalVariables;
      using Mon    = PhasorDynamics::Converter::ReecbMonitorableVariables;
      using Data   = PhasorDynamics::Converter::ReecbData<RealT, IdxT>;

      static constexpr size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(Ext variable)
      {
        return static_cast<size_t>(variable);
      }

      /// A model-vector row paired with an expected value. Rows are converted
      /// from `ReecbInternalVariables`, so the enum remains the single ordering
      /// contract shared with the implementation and README.
      using Row    = std::pair<size_t, RealT>;
      using Rows   = std::initializer_list<Row>;
      using ReecbT = PhasorDynamics::Converter::Reecb<ScalarT, IdxT>;

      /// A driven input value paired with the result it should produce. Kept
      /// distinct from `Row`, whose first member is a vector position.
      struct DrivenCase
      {
        RealT input;
        RealT expected;
      };

      /// Owns the terminal bus, REECB, assigned command nodes, and attached
      /// input nodes. Signal storage is declared before the model so every
      /// referenced node outlives REECB. Copying would invalidate the model and
      /// signal-node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, index(Ext::MAXIMUM)>                                   input_values_{};
        std::array<IdxT, index(Ext::MAXIMUM)>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, index(Ext::MAXIMUM)> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> iqcmd_node_;
        PhasorDynamics::SignalNode<T, IdxT> ipcmd_node_;

      public:
        explicit Fixture(const Data& data,
                         RealT       vr             = 1.0,
                         RealT       vi             = 0.0,
                         RealT       system_va_base = 100.0e6)
          : bus(static_cast<T>(vr), static_cast<T>(vi)),
            reecb(&bus, data)
        {
          reecb.setSystemBase(60.0, system_va_base);
          reecb.getSignals().template assignSignalNode<Vars::IQCMD>(&iqcmd_node_);
          reecb.getSignals().template assignSignalNode<Vars::IPCMD>(&ipcmd_node_);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        /// Attach fixture-owned storage to every external input.
        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = reecb.size() + bus.size();

          for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
          {
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = reecb.getSignals();
          signals.template attachSignalNode<Ext::PE>(&input_nodes_[index(Ext::PE)]);
          signals.template attachSignalNode<Ext::QGEN>(&input_nodes_[index(Ext::QGEN)]);
          signals.template attachSignalNode<Ext::QEXT>(&input_nodes_[index(Ext::QEXT)]);
          signals.template attachSignalNode<Ext::PFAREF>(&input_nodes_[index(Ext::PFAREF)]);
          signals.template attachSignalNode<Ext::PREF>(&input_nodes_[index(Ext::PREF)]);
        }

        /// Seed the assigned command nodes on the system base.
        void seedCommands(RealT iqcmd, RealT ipcmd)
        {
          iqcmd_node_.init(static_cast<T>(iqcmd));
          ipcmd_node_.init(static_cast<T>(ipcmd));
        }

        /// Everything REECB initialization requires: allocation, verification,
        /// an initialized terminal bus, and initialized command nodes.
        bool prepare(RealT iqcmd, RealT ipcmd)
        {
          const bool success = (bus.allocate() == 0) && (reecb.allocate() == 0)
                               && (reecb.verify() == 0) && (bus.initialize() == 0);
          if (!success)
          {
            std::cout << "REECB fixture preparation failed\n";
            return false;
          }

          seedCommands(iqcmd, ipcmd);
          return true;
        }

        /// prepare() plus successful REECB initialization.
        bool initialize(RealT iqcmd, RealT ipcmd)
        {
          if (!prepare(iqcmd, ipcmd))
          {
            return false;
          }
          if (reecb.initialize() != 0)
          {
            std::cout << "REECB initialization failed\n";
            return false;
          }
          return true;
        }

        int evaluate()
        {
          return reecb.evaluateResidual();
        }

        T iqcmd() const
        {
          return iqcmd_node_.read();
        }

        T ipcmd() const
        {
          return ipcmd_node_.read();
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
        PhasorDynamics::Converter::Reecb<T, IdxT> reecb;
      };

      static constexpr RealT kStateVr = 0.9;
      static constexpr RealT kStateVi = 0.4;

      Data makeMinimalData() const
      {
        Data data;
        data.device_class          = "Reecb";
        data.disambiguation_string = "reecb_test";
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);
        return data;
      }

      Data makeExplicitDefaultData() const
      {
        auto data = makeMinimalData();

        // These are the documented defaults. Vref0 is the terminal-voltage
        // fallback for the probe bus used by defaultsMatchDocumentedValues().
        data.parameters[Params::mva]    = 100.0;
        data.parameters[Params::PfFlag] = false;
        data.parameters[Params::VFlag]  = false;
        data.parameters[Params::QFlag]  = false;
        data.parameters[Params::Pqflag] = false;
        data.parameters[Params::Trv]    = 0.02;
        data.parameters[Params::Tp]     = 0.0;
        data.parameters[Params::Vref0]  = 0.9848857801796105;
        data.parameters[Params::Vdip]   = 0.85;
        data.parameters[Params::Vup]    = 1.15;
        data.parameters[Params::dbd1]   = 0.0;
        data.parameters[Params::dbd2]   = 0.0;
        data.parameters[Params::kqv]    = 5.0;
        data.parameters[Params::Iql1]   = -1.1;
        data.parameters[Params::Iqh1]   = 1.1;
        data.parameters[Params::Qmax]   = 0.436;
        data.parameters[Params::Qmin]   = -0.436;
        data.parameters[Params::Kqp]    = 0.0;
        data.parameters[Params::Kqi]    = 0.1;
        data.parameters[Params::Vmax]   = 1.1;
        data.parameters[Params::Vmin]   = 0.9;
        data.parameters[Params::Kvp]    = 18.0;
        data.parameters[Params::Kvi]    = 5.0;
        data.parameters[Params::Tiq]    = 0.02;
        data.parameters[Params::Tpord]  = 0.02;
        data.parameters[Params::dPmax]  = 99.0;
        data.parameters[Params::dPmin]  = -99.0;
        data.parameters[Params::Pmax]   = 1.0;
        data.parameters[Params::Pmin]   = 0.0;
        data.parameters[Params::Imax]   = 1.3;
        return data;
      }

      Data makeData() const
      {
        auto data = makeMinimalData();

        data.parameters[Params::mva]    = 100.0;
        data.parameters[Params::PfFlag] = false;
        data.parameters[Params::VFlag]  = true;
        data.parameters[Params::QFlag]  = false;
        data.parameters[Params::Pqflag] = true;
        data.parameters[Params::Trv]    = 0.02;
        data.parameters[Params::Tp]     = 0.02;
        data.parameters[Params::Vref0]  = 1.0;
        data.parameters[Params::Vdip]   = 0.7;
        data.parameters[Params::Vup]    = 1.2;
        data.parameters[Params::dbd1]   = -0.01;
        data.parameters[Params::dbd2]   = 0.01;
        data.parameters[Params::kqv]    = 0.0;
        data.parameters[Params::Iql1]   = -1.0;
        data.parameters[Params::Iqh1]   = 1.0;
        data.parameters[Params::Qmax]   = 1.0;
        data.parameters[Params::Qmin]   = -1.0;
        data.parameters[Params::Kqp]    = 1.0;
        data.parameters[Params::Kqi]    = 0.0;
        data.parameters[Params::Vmax]   = 1.2;
        data.parameters[Params::Vmin]   = 0.8;
        data.parameters[Params::Kvp]    = 1.0;
        data.parameters[Params::Kvi]    = 0.0;
        data.parameters[Params::Tiq]    = 0.02;
        data.parameters[Params::Tpord]  = 0.02;
        data.parameters[Params::dPmax]  = 1.0;
        data.parameters[Params::dPmin]  = -1.0;
        data.parameters[Params::Pmax]   = 1.0;
        data.parameters[Params::Pmin]   = 0.0;
        data.parameters[Params::Imax]   = 2.0;
        return data;
      }

      Data makeDynamicData() const
      {
        auto data                       = makeData();
        data.parameters[Params::mva]    = 50.0;
        data.parameters[Params::PfFlag] = true;
        data.parameters[Params::VFlag]  = true;
        data.parameters[Params::QFlag]  = true;
        data.parameters[Params::Pqflag] = true;
        data.parameters[Params::Trv]    = 0.2;
        data.parameters[Params::Tp]     = 0.4;
        data.parameters[Params::Vref0]  = 1.02;
        data.parameters[Params::Vdip]   = 0.7;
        data.parameters[Params::Vup]    = 1.2;
        data.parameters[Params::dbd1]   = -0.02;
        data.parameters[Params::dbd2]   = 0.03;
        data.parameters[Params::kqv]    = 2.0;
        data.parameters[Params::Iql1]   = -0.4;
        data.parameters[Params::Iqh1]   = 0.5;
        data.parameters[Params::Qmax]   = 0.8;
        data.parameters[Params::Qmin]   = -0.7;
        data.parameters[Params::Kqp]    = 0.6;
        data.parameters[Params::Kqi]    = 0.4;
        data.parameters[Params::Vmax]   = 1.3;
        data.parameters[Params::Vmin]   = 0.7;
        data.parameters[Params::Kvp]    = 1.2;
        data.parameters[Params::Kvi]    = 0.5;
        data.parameters[Params::Tiq]    = 0.3;
        data.parameters[Params::Tpord]  = 0.25;
        data.parameters[Params::dPmax]  = 0.6;
        data.parameters[Params::dPmin]  = -0.5;
        data.parameters[Params::Pmax]   = 1.4;
        data.parameters[Params::Pmin]   = 0.1;
        data.parameters[Params::Imax]   = 1.5;
        return data;
      }

      /// The external inputs the residual answer key is evaluated against.
      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::PE)     = 0.3;
        fixture.input(Ext::QGEN)   = -0.1;
        fixture.input(Ext::QEXT)   = 0.2;
        fixture.input(Ext::PFAREF) = 0.15;
        fixture.input(Ext::PREF)   = 0.35;
      }

      /// The rich state shared by the residual answer key and the Jacobian
      /// comparison. Every row is distinct so a swapped index cannot pass.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Converter::Reecb<T, IdxT>& reecb) const
      {
        setState(reecb,
                 {{index(Vars::VMEAS), 0.95},
                  {index(Vars::PMEAS), 0.55},
                  {index(Vars::XPIQ), 0.10},
                  {index(Vars::XPIV), -0.05},
                  {index(Vars::QV), 0.30},
                  {index(Vars::PORD), 0.65},
                  {index(Vars::VT), 1.00},
                  {index(Vars::VMEASSAFE), 0.96},
                  {index(Vars::SDIP), 0.80},
                  {index(Vars::VERR), 0.04},
                  {index(Vars::IQV), 0.15},
                  {index(Vars::QREF), 0.35},
                  {index(Vars::EQ), 0.20},
                  {index(Vars::VPIQ), 0.80},
                  {index(Vars::EPIV), -0.10},
                  {index(Vars::FPORD), 0.30},
                  {index(Vars::RPORD), 0.25},
                  {index(Vars::IQCIRC), 1.10},
                  {index(Vars::IPCIRC), 1.20},
                  {index(Vars::IQMAX), 1.00},
                  {index(Vars::IPMAX), 1.30},
                  {index(Vars::IQBASE), 0.20},
                  {index(Vars::IQRAW), 0.40},
                  {index(Vars::IQCMD), 0.25},
                  {index(Vars::IPCMD), 0.40}});
        setDerivative(reecb,
                      {{index(Vars::VMEAS), 0.01},
                       {index(Vars::PMEAS), -0.02},
                       {index(Vars::XPIQ), 0.03},
                       {index(Vars::XPIV), -0.04},
                       {index(Vars::QV), 0.05},
                       {index(Vars::PORD), -0.06}});
      }

      /// Omitting every optional parameter must give exactly the model built
      /// from the defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        auto implicit_data                    = makeMinimalData();
        implicit_data.parameters[Params::mva] = 100.0;

        Fixture<ScalarT> implicit_defaults(implicit_data, 0.9, 0.4);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), 0.9, 0.4);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        bool success = implicit_defaults.initialize(0.1, 0.2)
                       && explicit_defaults.initialize(0.1, 0.2);
        if (!success)
        {
          std::cout << "REECB documented-default comparison failed to initialize\n";
          return false;
        }

        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.reecb.y(),
                                   copyVector(explicit_defaults.reecb.y()),
                                   "documented-default state");
        success *= vectorUnchanged(implicit_defaults.reecb.yp(),
                                   copyVector(explicit_defaults.reecb.yp()),
                                   "documented-default derivative");
        success *= vectorUnchanged(implicit_defaults.reecb.getResidual(),
                                   copyVector(explicit_defaults.reecb.getResidual()),
                                   "documented-default residual");

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.reecb);
        setAnswerKeyState(explicit_defaults.reecb);
        success *= (implicit_defaults.evaluate() == 0);
        success *= (explicit_defaults.evaluate() == 0);
        success *= vectorUnchanged(implicit_defaults.reecb.getResidual(),
                                   copyVector(explicit_defaults.reecb.getResidual()),
                                   "documented-default dynamic residual");
        return success;
      }

      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus,
                                Params                              parameter,
                                RealT                               value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> model(&bus, data);
        return model.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected(PhasorDynamics::Bus<ScalarT, IdxT>& bus) const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>       unlinked_node;
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> model(&bus, makeData());
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

      bool initializationRejectedAtomically(const Data& data,
                                            RealT       iqcmd,
                                            RealT       ipcmd,
                                            RealT       terminal_voltage,
                                            const char* label) const
      {
        Fixture<ScalarT> fixture(data, terminal_voltage);
        fixture.attachAllInputs(77.0);
        if (!fixture.prepare(iqcmd, ipcmd))
        {
          return false;
        }

        auto* y  = fixture.reecb.y().getData();
        auto* yp = fixture.reecb.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(fixture.reecb.y().getSize()); ++i)
        {
          y[i]  = 0.125 + 0.01 * static_cast<RealT>(i);
          yp[i] = -0.25 - 0.01 * static_cast<RealT>(i);
        }
        fixture.seedCommands(iqcmd, ipcmd);
        fixture.reecb.y().setDataUpdated();
        fixture.reecb.yp().setDataUpdated();

        const auto y_before  = copyVector(fixture.reecb.y());
        const auto yp_before = copyVector(fixture.reecb.yp());

        bool success = true;
        if (fixture.reecb.initialize() == 0)
        {
          std::cout << "Expected initialization rejection: " << label << "\n";
          success = false;
        }

        success *= scalarMatches(fixture.iqcmd(), iqcmd, "rejected iqcmd preservation");
        success *= scalarMatches(fixture.ipcmd(), ipcmd, "rejected ipcmd preservation");
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
          success &= rowMatches(fixture.input(static_cast<Ext>(port)),
                                77.0,
                                "external input",
                                port,
                                "changed");
        }
        success *= vectorUnchanged(fixture.reecb.y(), y_before, "state");
        success *= vectorUnchanged(fixture.reecb.yp(), yp_before, "derivative");
        return success;
      }

      /// Write state rows and publish the update, folding in the
      /// setDataUpdated() that a hand-written write block has to remember.
      template <typename T>
      void setState(PhasorDynamics::Converter::Reecb<T, IdxT>& reecb, Rows rows) const
      {
        auto* y = reecb.y().getData();
        for (const auto& [row, value] : rows)
        {
          y[row] = static_cast<T>(value);
        }
        reecb.y().setDataUpdated();
      }

      /// setState() for the derivative vector.
      template <typename T>
      void setDerivative(PhasorDynamics::Converter::Reecb<T, IdxT>& reecb, Rows rows) const
      {
        auto* yp = reecb.yp().getData();
        for (const auto& [row, value] : rows)
        {
          yp[row] = static_cast<T>(value);
        }
        reecb.yp().setDataUpdated();
      }

      /// Compare one vector row against its expected value. Every row check in
      /// this suite reports through here, so failures share one format. Rows
      /// are named by the `ReecbInternalVariables` position the expectation was
      /// written with, leaving no parallel name string to maintain.
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
        std::cout << "REECB " << what << " row " << row << ' ' << context
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
                     const char*    context,
                     RealT          tolerance = kBehaviorTol) const
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
                                context,
                                tolerance);
        }
        return success;
      }

      bool residualsMatch(const ReecbT& reecb,
                          Rows          rows,
                          const char*   context   = "",
                          RealT         tolerance = kBehaviorTol) const
      {
        return rowsMatch(reecb.getResidual(),
                         rows.begin(),
                         rows.size(),
                         "residual",
                         context,
                         tolerance);
      }

      template <size_t size>
      bool residualsMatch(const ReecbT&                reecb,
                          const std::array<Row, size>& rows,
                          const char*                  context   = "",
                          RealT                        tolerance = kBehaviorTol) const
      {
        return rowsMatch(reecb.getResidual(), rows.data(), size, "residual", context, tolerance);
      }

      bool stateMatches(const ReecbT& reecb, Rows rows, const char* context = "") const
      {
        return rowsMatch(reecb.y(), rows.begin(), rows.size(), "state", context);
      }

      /// The model sits at a steady state: every residual and every derivative
      /// is zero.
      bool allResidualsZero(const ReecbT& reecb) const
      {
        bool        success = true;
        const auto* f       = reecb.getResidual().getData();
        const auto* yp      = reecb.yp().getData();
        for (size_t row = 0; row < static_cast<size_t>(reecb.getResidual().getSize()); ++row)
        {
          success &= rowMatches(static_cast<RealT>(f[row]),
                                0.0,
                                "residual",
                                row,
                                "at rest",
                                kSteadyStateTol);
          success &= rowMatches(static_cast<RealT>(yp[row]),
                                0.0,
                                "derivative",
                                row,
                                "at rest",
                                kSteadyStateTol);
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
      using DependencyMap = DependencyTracking::Variable::DependencyMap;

      static constexpr size_t externalVariableIndex(Ext variable)
      {
        constexpr size_t terminal_bus_size = 2;
        return index(Vars::MAXIMUM) + terminal_bus_size + index(variable);
      }

      bool dependencyMatches(const std::vector<DependencyMap>& jacobian,
                             Vars                              row,
                             size_t                            column,
                             RealT                             expected,
                             const char*                       context) const
      {
        const auto& dependencies = jacobian[index(row)];
        const auto  entry        = dependencies.find(column);
        RealT       actual       = 0.0;
        if (entry != dependencies.end())
        {
          actual = entry->second;
        }
        if (isEqual(actual, expected, kJacobianTol))
        {
          return true;
        }

        std::cout << "REECB Jacobian derivative (" << index(row) << ", " << column
                  << ") mismatch for " << context << ": " << std::setprecision(16)
                  << actual << " != " << expected << '\n';
        return false;
      }

      bool dependencyMatches(const std::vector<DependencyMap>& jacobian,
                             Vars                              row,
                             Vars                              column,
                             RealT                             expected,
                             const char*                       context) const
      {
        return dependencyMatches(jacobian, row, index(column), expected, context);
      }

      bool dependencyMatches(const std::vector<DependencyMap>& jacobian,
                             Vars                              row,
                             Ext                               column,
                             RealT                             expected,
                             const char*                       context) const
      {
        return dependencyMatches(jacobian,
                                 row,
                                 externalVariableIndex(column),
                                 expected,
                                 context);
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture) const
      {
        auto* y     = fixture.reecb.y().getData();
        auto* yp    = fixture.reecb.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        const auto model_size = static_cast<size_t>(fixture.reecb.size());
        for (size_t i = 0; i < model_size; ++i)
        {
          y[i].setVariableNumber(i);
          yp[i].setVariableNumber(i);
        }
        for (size_t i = 0; i < static_cast<size_t>(fixture.bus.size()); ++i)
        {
          bus_y[i].setVariableNumber(model_size + i);
        }
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
          const auto variable = static_cast<Ext>(port);
          fixture.input(variable).setVariableNumber(fixture.inputIndex(variable));
        }

        fixture.reecb.y().setDataUpdated();
        fixture.reecb.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          const Data& data,
          TestStatus& success) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.initialize(0.1, 0.2);
        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.reecb);
        numberVariables(fixture);
        success *= (fixture.evaluate() == 0);

        const auto                         model_size = static_cast<size_t>(fixture.reecb.size());
        std::vector<DepVar::DependencyMap> rows(model_size);
        const auto*                        f = fixture.reecb.getResidual().getData();
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
        success *= fixture.initialize(0.1, 0.2);

        for (IdxT i = 0; i < fixture.bus.size(); ++i)
        {
          fixture.bus.setVariableIndex(i, fixture.reecb.size() + i);
        }

        setAnswerKeyInputs(fixture);
        setAnswerKeyState(fixture.reecb);
        fixture.reecb.updateTime(0.0, 1.0);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.reecb.evaluateJacobian() == 0);
        success *= (fixture.reecb.constructCsr() == 0);
        return MapFromCsr(fixture.reecb.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
