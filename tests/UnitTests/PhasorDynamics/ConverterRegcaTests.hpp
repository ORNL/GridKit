#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    class ConverterRegcaTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterRegcaTests()  = default;
      ~ConverterRegcaTests() = default;

      static constexpr ScalarT kTol =
          static_cast<ScalarT>(100.0) * std::numeric_limits<ScalarT>::epsilon();

      /// Construction, the monitor, and every verify() error class: missing
      /// and invalid parameters, a null bus, and an unlinked command port.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> minimal(&bus);
        success *= (minimal.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (minimal.getMonitor() == nullptr);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> configured(&bus, makeData());
        success *= (configured.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() == 0);

        noteExpectedLogs("Testing invalid REGCA configurations. "
                         "Logged errors and time-constant warnings are expected.");
        success *= (minimal.verify() > 0);

        for (const Params parameter : {Params::Tg, Params::p0, Params::q0})
        {
          auto data = makeData();
          data.parameters.erase(parameter);
          PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing(&bus, data);
          success *= (missing.verify() > 0);
        }

        auto bad_switch                   = makeData();
        bad_switch.parameters[Params::sL] = static_cast<IdxT>(2);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> bad_switch_model(&bus, bad_switch);
        success *= (bad_switch_model.verify() > 0);

        success *= invalidParameterCase(bus, Params::mva, 0.0);
        success *= invalidParameterCase(bus, Params::IL1, -0.1);
        success *= invalidParameterCase(bus, Params::VL1, 0.3);
        success *= invalidParameterCase(bus, Params::VA1, 0.3);

        // Vhvmax must lie strictly above VA1: rejected at VA1 and just below.
        success *= invalidParameterCase(bus, Params::Vhvmax, kVa1);
        success *= invalidParameterCase(bus, Params::Vhvmax, kJustBelowVa1);

        // A null bus and an attached command with no linked source count as
        // configuration errors on the same footing as bad parameters.
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> busless(nullptr, makeData());
        success *= (busless.verify() > 0);

        PhasorDynamics::SignalNode<ScalarT, IdxT>       unlinked_node;
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> unlinked(&bus, makeData());
        unlinked.getSignals().template attachSignalNode<Ext::IPCMD>(&unlinked_node);
        success *= (unlinked.verify() > 0);

        // Zero time constants are raised to the well-posedness floor with a
        // warning, and the raised model still initializes to zero residuals.
        auto zero_time                   = makeData();
        zero_time.parameters[Params::Tg] = 0.0;
        zero_time.parameters[Params::TM] = 0.0;

        Fixture<ScalarT> fixture(zero_time);
        success *= fixture.initialize();
        success *= (fixture.evaluate() == 0);
        success *= allResidualsZero(fixture.regca);

        return success.report(__func__);
      }

      /// Power-flow initialization with every port exercised, then the
      /// fallback contract: unattached commands latch the initialization
      /// values as constant setpoints.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        // A 50 MVA converter on the 100 MVA system base at a (0.8, 0.6)
        // terminal. The command values are sentinels; initialize() must
        // publish the resolved steady-state commands over them.
        auto data                    = makeData();
        data.parameters[Params::mva] = 50.0;
        data.parameters[Params::p0]  = 0.4;
        data.parameters[Params::q0]  = -0.1;

        Fixture<ScalarT> fixture(data, 0.8, 0.6);
        fixture.attachIpcmd(-1.0);
        fixture.attachIqcmd(1.0);

        // Outputs alias y directly, so one published port pins the wiring.
        PhasorDynamics::SignalNode<ScalarT, IdxT> pbranch_node;
        fixture.regca.getSignals().template assignSignalNode<Vars::PBR>(&pbranch_node);

        success *= fixture.initialize();
        success *= (fixture.evaluate() == 0);

        // Fixed answer key for Vr = 0.8, Vi = 0.6, P0 = 0.4, Q0 = -0.1, and
        // a 50 MVA component base: IP = 2 P0 / VT and IQ = 2 Q0 / VT on the
        // component base, with IR and II returned to the system base.
        const auto* y  = fixture.regca.y().getData();
        success       *= scalarMatches(y[index(Vars::VM)], 1.0, "VM");
        success       *= scalarMatches(y[index(Vars::VT)], 1.0, "VT");
        success       *= scalarMatches(y[index(Vars::IP)], 0.80000000000025173, "IP");
        success       *= scalarMatches(y[index(Vars::IQ)], -0.2, "IQ");
        success       *= scalarMatches(y[index(Vars::IR)], 0.26, "IR");
        success       *= scalarMatches(y[index(Vars::II)], 0.32, "II");
        success       *= scalarMatches(y[index(Vars::PBR)], 0.4, "PBR");
        success       *= scalarMatches(y[index(Vars::QBR)], -0.1, "QBR");

        success *= scalarMatches(fixture.ipcmd, 0.40000000000012587, "published ipcmd");
        success *= scalarMatches(fixture.iqcmd, -0.1, "published iqcmd");
        success *= scalarMatches(pbranch_node.read(), 0.4, "pbranch signal");

        // The accumulated bus injection is the system-base branch current.
        success *= scalarMatches(fixture.bus.Ir(), 0.26, "bus real injection");
        success *= scalarMatches(fixture.bus.Ii(), 0.32, "bus imaginary injection");

        success *= allResidualsZero(fixture.regca);

        // Without attached ports the commands latch p0 and q0. Displace the
        // current states by different amounts; the latched commands pull
        // both back at the interior first-order rates.
        auto latch_data                   = makeDynamicData();
        latch_data.parameters[Params::q0] = 0.2;

        Fixture<ScalarT> latched(latch_data);
        success *= latched.initialize();

        auto* latched_y  = latched.regca.y().getData();
        success         *= scalarMatches(latched_y[index(Vars::IP)],
                                 0.60000000000018872,
                                 "initialized IP");
        success         *= scalarMatches(latched_y[index(Vars::IQ)], 0.2, "initialized IQ");

        latched_y[index(Vars::IP)] = 0.5;
        latched_y[index(Vars::IQ)] = 0.14;
        latched.regca.y().setDataUpdated();
        success *= (latched.evaluate() == 0);

        // f[IP] = (0.6 - 0.5) / Tg and f[IQ] = (0.2 - 0.14) / Tg with
        // Tg = 0.2; both rates sit inside every limiter.
        const auto* f  = latched.regca.getResidual().getData();
        success       *= scalarMatches(f[index(Vars::IP)],
                                 0.50000000000094358,
                                 "latched active-current rate");
        success       *= scalarMatches(f[index(Vars::IQ)], 0.3, "latched reactive-current rate");

        return success.report(__func__);
      }

      /// The admissible initialization domain: every rejected operating
      /// point, then the accepted boundaries next to them.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible REGCA initialization points. "
                         "Logged errors are expected.");

        struct RejectionCase
        {
          const char* label;
          RealT       vr;
          RealT       p0;
          RealT       il1;
        };

        // P0 and IL1 default to the makeData() values; each row breaks
        // exactly one initialize() guard.
        const std::array<RejectionCase, 4> rejected{{
            {"terminal voltage at Vhvmax", kHvrcmVoltageLimit, 0.0, 1.1},
            {"terminal voltage above Vhvmax", kHvrcmVoltageLimit + 0.1, 0.0, 1.1},
            {"terminal voltage just below VA1", kJustBelowVa1, 0.0, 1.1},
            {"zero terminal voltage", 0.0, 0.0, 1.1},
        }};

        for (const auto& test_case : rejected)
        {
          auto data                    = makeData();
          data.parameters[Params::p0]  = test_case.p0;
          data.parameters[Params::IL1] = test_case.il1;

          Fixture<ScalarT> fixture(data, test_case.vr);
          success *= fixture.prepare();
          if (fixture.regca.initialize() == 0)
          {
            std::cout << "Expected rejection: " << test_case.label << "\n";
            success = false;
          }
        }

        // A finite ceiling binds only inside the LVPL ramp: with VL1 raised
        // above the terminal voltage, an active current beyond the ramp value
        // is rejected.
        {
          auto data                    = makeData();
          data.parameters[Params::p0]  = 0.6;
          data.parameters[Params::IL1] = 0.2;
          data.parameters[Params::VL1] = 1.05;

          Fixture<ScalarT> fixture(data);
          success *= fixture.prepare();
          if (fixture.regca.initialize() == 0)
          {
            std::cout << "Expected rejection: LVPL ceiling below the active current\n";
            success = false;
          }
        }

        // VA1 itself is an admissible LVACM boundary and the operating
        // point reproduces the P0/Q0 injection.
        {
          auto data                   = makeData();
          data.parameters[Params::p0] = 0.2;
          data.parameters[Params::q0] = 0.1;

          Fixture<ScalarT> fixture(data, kVa1);
          success *= fixture.initialize();
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.regca);

          const auto* y  = fixture.regca.y().getData();
          success       *= scalarMatches(y[index(Vars::PBR)], 0.2, "PBR at the LVACM upper breakpoint");
          success       *= scalarMatches(y[index(Vars::QBR)], 0.1, "QBR at the LVACM upper breakpoint");
        }

        // A ceiling just above the requested current is admissible.
        {
          auto data                    = makeData();
          data.parameters[Params::p0]  = 0.6;
          data.parameters[Params::IL1] = 0.61;

          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize();
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.regca);

          const auto* y  = fixture.regca.y().getData();
          success       *= scalarMatches(y[index(Vars::IP)],
                                   0.60000000000018872,
                                   "IP below the LVPL ceiling");
          success       *= scalarMatches(y[index(Vars::PBR)], 0.6, "PBR below the LVPL ceiling");
        }

        // The active-current state may initialize exactly on the LVPL ceiling.
        // VA1 is lowered so the terminal voltage sits inside the LVPL ramp,
        // where the release term vanishes.
        {
          auto data                    = makeData();
          data.parameters[Params::IL1] = 0.0;
          data.parameters[Params::VA1] = 0.5;

          Fixture<ScalarT> fixture(data, 0.6);
          success *= fixture.initialize();
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.regca);

          const auto* y  = fixture.regca.y().getData();
          success       *= scalarMatches(y[index(Vars::IP)], 0.0, "IP at the LVPL ceiling");
          success       *= scalarMatches(y[index(Vars::IL)], 0.0, "LVPL ceiling at IP");
        }

        // Above the upper breakpoint the ceiling releases: an operating point
        // beyond IL1 initializes at healthy voltage.
        {
          auto data                   = makeData();
          data.parameters[Params::p0] = 1.3;

          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize();
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.regca);

          const auto* y  = fixture.regca.y().getData();
          success       *= scalarMatches(y[index(Vars::IP)],
                                   1.3000000000004091,
                                   "IP beyond IL1 with the ceiling released");
          success       *= scalarMatches(y[index(Vars::PBR)], 1.3, "PBR beyond IL1 with the ceiling released");
        }

        // With LVPL bypassed a collapsed ceiling does not constrain the
        // initialization.
        {
          auto data                    = makeData();
          data.parameters[Params::p0]  = 0.6;
          data.parameters[Params::sL]  = false;
          data.parameters[Params::IL1] = 0.2;

          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize();
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.regca);

          const auto* y  = fixture.regca.y().getData();
          success       *= scalarMatches(y[index(Vars::IP)],
                                   0.60000000000018872,
                                   "IP with LVPL bypassed");
        }

        return success.report(__func__);
      }

      /// The complete equation answer key at a hand-computed state.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeDynamicData(), kStateVr, kStateVi);
        fixture.attachIpcmd(kStateIpcmd);
        fixture.attachIqcmd(kStateIqcmd);
        success *= fixture.initialize();

        // initialize() published steady-state commands over the attached
        // values; restore the commands the answer key assumes.
        fixture.ipcmd = kStateIpcmd;
        fixture.iqcmd = kStateIqcmd;

        setAnswerKeyState(fixture.regca);
        success *= (fixture.evaluate() == 0);

        struct ExpectedResidual
        {
          Vars        variable;
          const char* name;
          RealT       value;
        };

        // Each entry is the expected value of the named README equation at
        // the answer-key state.
        const std::array<ExpectedResidual, 10> expected{{
            {Vars::VM, "VM", 0.865},               // -VM' + (VT - VM) / TM
            {Vars::IQ, "IQ", 1.52},                // -IQ' + max(fq, Rqmin)
            {Vars::IP, "IP", 0.22},                // -IP' + IL' + awmax(IP - IL, fp_limited - IL', 0)
            {Vars::VT, "VT", -0.035},              // -VT^2 + Vr^2 + Vi^2
            {Vars::IR, "IR", 0.24999999999974593}, // -VT*IR + Vi*(IQ - IQEXTRA) + Vr*IP*linseg(VT)
            {Vars::II, "II", 0.25099999999993317}, // -VT*II - Vr*(IQ - IQEXTRA) + Vi*IP*linseg(VT)
            {Vars::IQEXTRA, "IQEXTRA", -0.03},     // smooth HVRCM constraint
            {Vars::IL, "IL", 0.35},                // -IL + linseg(VM, VL0, VL1, IL1)
            {Vars::PBR, "PBR", 0.0},               // -PBR + Vr*IR + Vi*II
            {Vars::QBR, "QBR", 0.0},               // -QBR + Vi*IR - Vr*II
        }};

        const auto& residual  = fixture.regca.getResidual();
        success              *= (static_cast<size_t>(residual.getSize()) == expected.size());

        const auto* f = residual.getData();
        for (const auto& row : expected)
        {
          success *= scalarMatches(f[index(row.variable)], row.value, row.name);
        }

        return success.report(__func__);
      }

      /// Active-current control: rrpwr limits motion that increases |Ip|,
      /// smoothly releases restoring motion, and LVPL bounds the integrator
      /// only when enabled.
      TestOutcome activeCurrentControl()
      {
        TestStatus success = true;

        // With LVPL bypassed, select commands that drive Ip outward or toward
        // zero. At Ip = 0, rrpwr applies the symmetric Rpmax = 0.7 limit in
        // either direction.
        {
          auto data                   = makeDynamicData();
          data.parameters[Params::sL] = false;

          Fixture<ScalarT> fixture(data);
          fixture.attachIpcmd(0.0);
          success *= fixture.initialize();

          struct RateCase
          {
            const char* label;
            RealT       current;
            RealT       command;
            RealT       expected_rate;
          };

          const std::array<RateCase, 6> cases{{
              {"negative outward rate", -0.5, -1.5, -0.7},
              {"negative restoring rate", -0.5, 0.5, 5.0},
              {"zero negative rate", 0.0, -0.5, -0.7},
              {"zero positive rate", 0.0, 0.5, 0.7},
              {"positive restoring rate", 0.5, -0.5, -5.0},
              {"positive outward rate", 0.5, 1.5, 0.7},
          }};

          auto* y = fixture.regca.y().getData();
          for (const auto& test_case : cases)
          {
            fixture.ipcmd      = test_case.command;
            y[index(Vars::IP)] = test_case.current;
            fixture.regca.y().setDataUpdated();

            success *= (fixture.evaluate() == 0);

            const auto* f  = fixture.regca.getResidual().getData();
            success       *= scalarMatches(f[index(Vars::IP)],
                                     test_case.expected_rate,
                                     test_case.label);
          }
        }

        struct LvplCase
        {
          const char* label;
          bool        enabled;
          RealT       command;
          RealT       expected_rate;
        };

        // Above the LVPL ceiling, outward motion is blocked and restoring
        // motion passes. Bypassing LVPL leaves the same outward drive intact.
        const std::array<LvplCase, 3> cases{{
            {"outward rate blocked by LVPL", true, 0.7, 0.0},
            {"restoring rate passes LVPL", true, 0.5, -0.5},
            {"outward rate with LVPL bypassed", false, 0.7, 0.5},
        }};

        for (const auto& test_case : cases)
        {
          auto data                   = makeDynamicData();
          data.parameters[Params::sL] = test_case.enabled;

          Fixture<ScalarT> fixture(data);
          fixture.attachIpcmd(test_case.command);
          success *= fixture.initialize();

          // initialize() publishes the steady-state command; restore the
          // driven value before evaluating the limiter behavior.
          fixture.ipcmd = test_case.command;

          auto* y            = fixture.regca.y().getData();
          y[index(Vars::IP)] = 0.6;
          y[index(Vars::IL)] = 0.4;
          fixture.regca.y().setDataUpdated();

          success *= (fixture.evaluate() == 0);

          const auto* f  = fixture.regca.getResidual().getData();
          success       *= scalarMatches(f[index(Vars::IP)],
                                   test_case.expected_rate,
                                   test_case.label);
        }

        // The moving-ceiling rule: with VM inside the active LVPL segment and
        // the sensed voltage falling, a pinned Ip tracks the ceiling rate
        // downward instead of freezing. IL' = 2.2 * (0.3 - 0.65) / 0.4.
        {
          Fixture<ScalarT> fixture(makeDynamicData());
          fixture.attachIpcmd(0.7);
          success *= fixture.initialize();

          fixture.ipcmd = 0.7;

          auto* y            = fixture.regca.y().getData();
          y[index(Vars::IP)] = 0.6;
          y[index(Vars::IL)] = 0.4;
          y[index(Vars::VM)] = 0.65;
          y[index(Vars::VT)] = 0.3;
          fixture.regca.y().setDataUpdated();

          success *= (fixture.evaluate() == 0);

          const auto* f  = fixture.regca.getResidual().getData();
          success       *= scalarMatches(f[index(Vars::IP)],
                                   -1.925,
                                   "falling ceiling drags Ip");
        }

        // Above the upper breakpoint the ceiling releases: outward motion
        // beyond IL1 passes at healthy voltage, limited only by rrpwr.
        {
          Fixture<ScalarT> fixture(makeDynamicData());
          fixture.attachIpcmd(1.5);
          success *= fixture.initialize();

          fixture.ipcmd = 1.5;

          auto* y            = fixture.regca.y().getData();
          y[index(Vars::IP)] = 1.3;
          fixture.regca.y().setDataUpdated();

          success *= (fixture.evaluate() == 0);

          const auto* f  = fixture.regca.getResidual().getData();
          success       *= scalarMatches(f[index(Vars::IP)],
                                   0.7,
                                   "released ceiling above the breakpoint");
        }

        return success.report(__func__);
      }

      /// The recovery-rate limiter branch is chosen by the sign of Q0 at
      /// initialization, not by the live command: positive Q0 caps rising
      /// rates at Rqmax, negative Q0 floors falling rates at Rqmin, and zero
      /// Q0 leaves both directions unrestricted.
      TestOutcome reactiveCurrentControl()
      {
        TestStatus success = true;

        struct RecoveryCase
        {
          const char* label;
          RealT       q0;
          RealT       command;
          RealT       expected_rate;
        };

        // The unlimited rate is (command - Q0) / Tg with Tg = 0.2.
        const std::array<RecoveryCase, 4> cases{{
            {"upper recovery limit binds", 0.2, 0.4, 0.5},
            {"lower recovery limit binds", -0.2, -0.4, -0.6},
            {"zero Q0 leaves rising rates unrestricted", 0.0, 0.4, 2.0},
            {"zero Q0 leaves falling rates unrestricted", 0.0, -0.4, -2.0},
        }};

        for (const auto& test_case : cases)
        {
          auto data                   = makeDynamicData();
          data.parameters[Params::q0] = test_case.q0;

          Fixture<ScalarT> fixture(data);
          fixture.attachIqcmd(test_case.q0);
          success *= fixture.initialize();

          fixture.iqcmd  = test_case.command;
          success       *= (fixture.evaluate() == 0);

          const auto* f  = fixture.regca.getResidual().getData();
          success       *= scalarMatches(f[index(Vars::IQ)], test_case.expected_rate, test_case.label);
        }

        return success.report(__func__);
      }

      /// High-voltage reactive current management: the initialization root,
      /// the residual across the transition, and its local derivative.
      TestOutcome highVoltageManagement()
      {
        TestStatus success = true;

        // Initialization near Vhvmax solves the smooth constraint for a
        // nonzero IQEXTRA. Half the transition margin keeps the root
        // distinct from the margin, so a swapped-argument residual cannot
        // satisfy this operating point.
        {
          const RealT margin           = 0.5 * hvrcmTransition();
          const RealT terminal_voltage = kHvrcmVoltageLimit - margin;

          auto data                   = makeData();
          data.parameters[Params::q0] = 0.1;

          Fixture<ScalarT> fixture(data, terminal_voltage);
          success *= fixture.initialize();
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.regca);

          // The zero residual already pins the root to the HVRCM equation.
          // What remains to check is that it is a genuine compensation
          // above the margin it cancels, and that it leaves Q0 intact.
          const auto* y = fixture.regca.y().getData();
          if (!(y[index(Vars::IQEXTRA)] > margin))
          {
            std::cout << "IQEXTRA is not above the voltage margin it compensates\n";
            success = false;
          }
          success *= scalarMatches(y[index(Vars::IQ)] - y[index(Vars::IQEXTRA)],
                                   0.1 / terminal_voltage,
                                   "IQ preserves Q0 after HVRCM compensation");
          success *= scalarMatches(y[index(Vars::QBR)], 0.1, "QBR");
        }

        // Numeric answer keys pin both sign and magnitude of the HVRCM row.
        {
          Fixture<ScalarT> fixture(makeData());
          success *= fixture.initialize();

          struct HvrcmCase
          {
            RealT voltage;
            RealT extra_current;
            RealT expected_residual;
          };

          const RealT transition = hvrcmTransition();

          // Away from the transition the smooth tail is below the tolerance
          // and the row reduces to -IQEXTRA; at and beyond the ceiling the
          // tail is the whole answer.
          const std::array<HvrcmCase, 5> cases{{
              {kHvrcmVoltageLimit - 0.2, 0.05, -0.05},
              {kHvrcmVoltageLimit - 0.1, -0.05, 0.05},
              {kHvrcmVoltageLimit - transition, transition, 0.0},
              {kHvrcmVoltageLimit, 0.0, transition},
              {kHvrcmVoltageLimit + 0.1, 0.1, 0.1},
          }};

          auto* y = fixture.regca.y().getData();
          for (const auto& test_case : cases)
          {
            y[index(Vars::VT)]      = test_case.voltage;
            y[index(Vars::IQEXTRA)] = test_case.extra_current;
            fixture.regca.y().setDataUpdated();

            success *= (fixture.evaluate() == 0);

            const auto* f  = fixture.regca.getResidual().getData();
            success       *= scalarMatches(f[index(Vars::IQEXTRA)],
                                     test_case.expected_residual,
                                     "HVRCM residual");
          }
        }

        // At the transition point both partial derivatives of the HVRCM row
        // are exactly one half.
        {
          using DepVar = DependencyTracking::Variable;

          Fixture<DepVar> fixture(makeData());
          success *= fixture.initialize();

          auto* y                 = fixture.regca.y().getData();
          y[index(Vars::VT)]      = kHvrcmVoltageLimit - hvrcmTransition();
          y[index(Vars::IQEXTRA)] = hvrcmTransition();
          fixture.regca.y().setDataUpdated();
          numberVariables(fixture);

          success *= (fixture.evaluate() == 0);

          const auto& dependencies =
              fixture.regca.getResidual().getData()[index(Vars::IQEXTRA)].getDependencies();

          const DepVar::DependencyMap expected{{
              {index(Vars::VT), 0.5},
              {index(Vars::IQEXTRA), -0.5},
          }};
          success *= isEqual(dependencies, expected, kTol);
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Each LVPL configuration drives both sensitivity paths, so the
      /// compared rows cannot drift apart between hand-built fixtures.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const RealT                transition_current = 1.0 / Math::MU<RealT>;
        const std::array<RealT, 3> currents{{
            -transition_current,
            0.0,
            transition_current,
        }};

        for (const bool lvpl_enabled : {false, true})
        {
          auto data                   = makeJacobianData();
          data.parameters[Params::sL] = lvpl_enabled;

          for (const RealT current : currents)
          {
            const auto dependency_tracking_jacobian =
                dependencyTrackingJacobian(data, current, success);
            const auto enzyme_jacobian = enzymeJacobian(data, current, success);

            const auto ip_row  = index(Vars::IP);
            const auto il_col  = index(Vars::IL);
            success           *= dependency_tracking_jacobian[ip_row].contains(il_col);
            success           *= enzyme_jacobian[ip_row].contains(il_col);

            success          *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
            const auto nrows  = std::min(dependency_tracking_jacobian.size(),
                                        enzyme_jacobian.size());

            for (size_t i = 0; i < nrows; ++i)
            {
              if (!isEqual(dependency_tracking_jacobian[i],
                           enzyme_jacobian[i],
                           kTol))
              {
                std::cout << "Jacobian row " << i
                          << " mismatch between dependency tracking and Enzyme"
                          << " at IP = " << current << "\n";
                success = false;
              }
            }
          }
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Converter::RegcaParameters;
      using Vars   = PhasorDynamics::Converter::RegcaInternalVariables;
      using Ext    = PhasorDynamics::Converter::RegcaExternalVariables;
      using Data   = PhasorDynamics::Converter::RegcaData<RealT, IdxT>;

      /// Owns a terminal bus, the model under test, and its two command
      /// signals. Copying would dangle the bus pointer regca holds, so the
      /// fixture is pinned.
      template <typename T>
      struct Fixture
      {
        explicit Fixture(const Data& data, RealT vr = 1.0, RealT vi = 0.0)
          : bus(static_cast<T>(vr), static_cast<T>(vi)),
            regca(&bus, data),
            ipcmd_index(regca.size() + bus.size()),
            iqcmd_index(ipcmd_index + 1)
        {
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachIpcmd(RealT value)
        {
          ipcmd = value;
          ipcmd_node.set(&ipcmd, &ipcmd_index);
          regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        }

        void attachIqcmd(RealT value)
        {
          iqcmd = value;
          iqcmd_node.set(&iqcmd, &iqcmd_index);
          regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);
        }

        /// Everything initialize() requires: allocation, verification, and
        /// an initialized terminal bus.
        bool prepare()
        {
          const bool success = (bus.allocate() == 0) && (regca.allocate() == 0)
                               && (regca.verify() == 0) && (bus.initialize() == 0);
          if (!success)
          {
            std::cout << "REGCA fixture preparation failed\n";
          }
          return success;
        }

        /// prepare() plus a successful model initialization.
        bool initialize()
        {
          if (!prepare())
          {
            return false;
          }
          if (regca.initialize() != 0)
          {
            std::cout << "REGCA initialization failed\n";
            return false;
          }
          return true;
        }

        /// Zeroes the bus injection, then accumulates regca into it. The
        /// ordering the models require.
        int evaluate()
        {
          const int bus_status = bus.evaluateResidual();
          if (bus_status != 0)
          {
            return bus_status;
          }
          return regca.evaluateResidual();
        }

        PhasorDynamics::Bus<T, IdxT>              bus;
        PhasorDynamics::Converter::Regca<T, IdxT> regca;

        T    ipcmd{};
        T    iqcmd{};
        IdxT ipcmd_index;
        IdxT iqcmd_index;

        PhasorDynamics::SignalNode<T, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<T, IdxT> iqcmd_node;
      };

      static size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      // Bus and command values shared by the residual and Jacobian fixtures.
      static constexpr RealT kStateVr    = 0.95;
      static constexpr RealT kStateVi    = 0.25;
      static constexpr RealT kStateIpcmd = 0.9;
      static constexpr RealT kStateIqcmd = 0.1;

      static constexpr RealT kVa1          = 0.9;
      static constexpr RealT kJustBelowVa1 = kVa1 - 1.0e-6;

      /// Vhvmax in makeData() and in makeDynamicData().
      static constexpr RealT kHvrcmVoltageLimit        = 1.2;
      static constexpr RealT kDynamicHvrcmVoltageLimit = 1.3;

      /// The smooth ramp at a zero argument. Used as an HVRCM voltage
      /// margin, this is the one margin whose constraint root is the margin
      /// itself, so the HVRCM row is exactly zero there and its two partial
      /// derivatives are exactly one half.
      static RealT hvrcmTransition()
      {
        return Math::ramp(static_cast<RealT>(0.0));
      }

      Data makeData()
      {
        using Buses = PhasorDynamics::Converter::RegcaBuses;
        using Mon   = PhasorDynamics::Converter::RegcaMonitorableVariables;

        Data data;
        data.device_class          = "Regca";
        data.disambiguation_string = "regca_test";
        data.buses[Buses::bus]     = 1;
        data.monitored_variables.insert(Mon::ir);
        data.monitored_variables.insert(Mon::ii);
        data.monitored_variables.insert(Mon::p);
        data.monitored_variables.insert(Mon::q);

        data.parameters[Params::p0]     = 0.0;
        data.parameters[Params::q0]     = 0.0;
        data.parameters[Params::mva]    = 100.0;
        data.parameters[Params::Tg]     = 0.02;
        data.parameters[Params::TM]     = 0.02;
        data.parameters[Params::Rqmax]  = 999.0;
        data.parameters[Params::Rqmin]  = -999.0;
        data.parameters[Params::Rpmax]  = 999.0;
        data.parameters[Params::sL]     = true;
        data.parameters[Params::IL1]    = 1.1;
        data.parameters[Params::VL0]    = 0.4;
        data.parameters[Params::VL1]    = 0.9;
        data.parameters[Params::VA0]    = 0.4;
        data.parameters[Params::VA1]    = 0.9;
        data.parameters[Params::Vhvmax] = kHvrcmVoltageLimit;

        return data;
      }

      /// Dynamic-response parameters: slower lags and tight rate limits so
      /// limiter selections show up in the residual rows.
      Data makeDynamicData()
      {
        auto data                       = makeData();
        data.parameters[Params::p0]     = 0.6;
        data.parameters[Params::q0]     = -0.2;
        data.parameters[Params::Tg]     = 0.2;
        data.parameters[Params::TM]     = 0.4;
        data.parameters[Params::Rqmax]  = 0.5;
        data.parameters[Params::Rqmin]  = -0.6;
        data.parameters[Params::Rpmax]  = 0.7;
        data.parameters[Params::Vhvmax] = kDynamicHvrcmVoltageLimit;
        return data;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Jacobian data: dynamic parameters on a non-identity power base.
      /// P0 = 0 keeps the reduced-base active current below the LVPL
      /// ceiling at initialization.
      Data makeJacobianData()
      {
        auto data                    = makeDynamicData();
        data.parameters[Params::p0]  = 0.0;
        data.parameters[Params::mva] = 50.0;
        return data;
      }
#endif

      /// Writes the answer-key state shared by the residual and Jacobian
      /// tests. Values only; in the dependency-tracking fixture, number
      /// variables after this call.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Converter::Regca<T, IdxT>& regca)
      {
        auto* y                 = regca.y().getData();
        y[index(Vars::VM)]      = 0.65;
        y[index(Vars::IQ)]      = -0.2;
        y[index(Vars::IP)]      = 0.85;
        y[index(Vars::VT)]      = 1.0;
        y[index(Vars::IR)]      = 0.5;
        y[index(Vars::II)]      = 0.18;
        y[index(Vars::IQEXTRA)] = 0.03;
        y[index(Vars::IL)]      = 0.2;
        y[index(Vars::PBR)]     = 0.52;
        y[index(Vars::QBR)]     = -0.046;

        auto* yp            = regca.yp().getData();
        yp[index(Vars::VM)] = 0.01;
        yp[index(Vars::IQ)] = -0.02;
        yp[index(Vars::IP)] = 0.03;

        regca.y().setDataUpdated();
        regca.yp().setDataUpdated();
      }

      /// Numbers regca y and yp together, the bus after the regca block,
      /// and the command signals at their port indices, matching the
      /// Jacobian layout. Write state values first; numbering resets each
      /// dependency map. Numbering an unattached command is a no-op for the
      /// model.
      void numberVariables(Fixture<DependencyTracking::Variable>& fixture)
      {
        auto* y     = fixture.regca.y().getData();
        auto* yp    = fixture.regca.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        const auto regca_size = static_cast<size_t>(fixture.regca.size());
        for (size_t i = 0; i < regca_size; ++i)
        {
          y[i].setVariableNumber(i);
          yp[i].setVariableNumber(i);
        }
        for (size_t i = 0; i < static_cast<size_t>(fixture.bus.size()); ++i)
        {
          bus_y[i].setVariableNumber(i + regca_size);
        }
        fixture.ipcmd.setVariableNumber(fixture.ipcmd_index);
        fixture.iqcmd.setVariableNumber(fixture.iqcmd_index);

        fixture.regca.y().setDataUpdated();
        fixture.regca.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus, Params param, RealT value)
      {
        auto data              = makeData();
        data.parameters[param] = value;
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> model(&bus, data);
        return model.verify() > 0;
      }

      bool allResidualsZero(PhasorDynamics::Converter::Regca<ScalarT, IdxT>& regca) const
      {
        bool        success = true;
        const auto* f       = regca.getResidual().getData();
        const auto* yp      = regca.yp().getData();
        for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
        {
          if (!isEqual(f[i], 0.0, kTol))
          {
            std::cout << "REGCA residual row " << i << " is " << f[i] << "\n";
            success = false;
          }
          if (!isEqual(yp[i], 0.0, kTol))
          {
            std::cout << "REGCA derivative row " << i << " is " << yp[i] << "\n";
            success = false;
          }
        }
        return success;
      }

      bool scalarMatches(ScalarT actual, ScalarT expected, const char* label) const
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }

        std::cout << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << "\n";
        return false;
      }

      /// Raises verbosity just long enough to flag intentionally provoked
      /// logs, so expected errors in the output are not mistaken for
      /// failures.
      void noteExpectedLogs(const char* message) const
      {
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << "\n";
        Log::setVerbosity(previous_verbosity);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Move the answer-key state to the HVRCM and active-current limiter
      /// transition points, where the Jacobian has the richest structure.
      template <typename T>
      void setJacobianState(
          PhasorDynamics::Converter::Regca<T, IdxT>& regca,
          RealT                                      current)
      {
        setAnswerKeyState(regca);
        auto* y                 = regca.y().getData();
        y[index(Vars::IP)]      = current;
        y[index(Vars::IL)]      = 0.0;
        y[index(Vars::VT)]      = kDynamicHvrcmVoltageLimit - hvrcmTransition();
        y[index(Vars::IQEXTRA)] = hvrcmTransition();
        regca.y().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian(
          const Data& data,
          RealT       current,
          TestStatus& success)
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachIpcmd(kStateIpcmd);
        fixture.attachIqcmd(kStateIqcmd);
        success *= fixture.initialize();

        fixture.ipcmd = kStateIpcmd;
        fixture.iqcmd = kStateIqcmd;
        setJacobianState(fixture.regca, current);
        numberVariables(fixture);

        success *= (fixture.evaluate() == 0);

        const auto  regca_size = static_cast<size_t>(fixture.regca.size());
        const auto* f          = fixture.regca.getResidual().getData();

        std::vector<DepVar::DependencyMap> dependencies(
            regca_size + static_cast<size_t>(fixture.bus.size()));
        for (size_t i = 0; i < regca_size; ++i)
        {
          dependencies[i] = f[i].getDependencies();
        }
        dependencies[regca_size]     = fixture.bus.Ir().getDependencies();
        dependencies[regca_size + 1] = fixture.bus.Ii().getDependencies();

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian(
          const Data& data,
          RealT       current,
          TestStatus& success)
      {
        Fixture<ScalarT> fixture(data, kStateVr, kStateVi);
        fixture.attachIpcmd(kStateIpcmd);
        fixture.attachIqcmd(kStateIqcmd);
        success *= fixture.initialize();

        // Number the bus variables and residual rows after the regca block,
        // matching the dependency-tracking layout. evaluateJacobian() reads
        // these indices at call time, so assigning them after bring-up works.
        for (IdxT i = 0; i < fixture.bus.size(); ++i)
        {
          fixture.bus.setVariableIndex(i, i + fixture.regca.size());
          fixture.bus.setResidualIndex(i, i + fixture.regca.size());
        }

        fixture.ipcmd = kStateIpcmd;
        fixture.iqcmd = kStateIqcmd;
        setJacobianState(fixture.regca, current);
        fixture.regca.updateTime(0.0, 1.0);

        success *= (fixture.evaluate() == 0);
        success *= (fixture.regca.evaluateJacobian() == 0);
        success *= (fixture.regca.constructCsr() == 0);

        return MapFromCsr(fixture.regca.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
