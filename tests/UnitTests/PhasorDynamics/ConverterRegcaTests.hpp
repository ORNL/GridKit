#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
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
    class ConverterRegcaTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterRegcaTests()  = default;
      ~ConverterRegcaTests() = default;

      // At nominal voltage, the CommonMath LVACM approximation differs from
      // its piecewise plateau by O(1e-13).
      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-12);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> minimal(&bus);
        success *= (minimal.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (minimal.getMonitor() == nullptr);

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing incomplete REGCA configuration. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);
        success *= (minimal.verify() > 0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> configured(&bus, makeData());
        success *= (configured.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing invalid and missing REGCA parameters. "
                    << "Logged errors and time-constant warnings are expected.\n";
        Log::setVerbosity(previous_verbosity);

        auto missing = makeData();
        missing.parameters.erase(Params::Tg);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing_model(&bus, missing);
        success *= (missing_model.verify() > 0);

        auto missing_p0 = makeData();
        missing_p0.parameters.erase(Params::p0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing_p0_model(&bus, missing_p0);
        success *= (missing_p0_model.verify() > 0);

        auto missing_q0 = makeData();
        missing_q0.parameters.erase(Params::q0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing_q0_model(&bus, missing_q0);
        success *= (missing_q0_model.verify() > 0);

        auto bad_switch                   = makeData();
        bad_switch.parameters[Params::sL] = static_cast<IdxT>(2);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> bad_switch_model(&bus, bad_switch);
        success *= (bad_switch_model.verify() > 0);

        success *= invalidParameterCase(bus, Params::mva, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Rpmax, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Rqmin, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::IL1, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::VL1, static_cast<RealT>(0.3));
        success *= invalidParameterCase(bus, Params::VA1, static_cast<RealT>(0.3));
        success *= invalidParameterCase(bus, Params::Vhvmax, static_cast<RealT>(0.0));

        auto zero_time_constants                   = makeData();
        zero_time_constants.parameters[Params::Tg] = static_cast<RealT>(0.0);
        zero_time_constants.parameters[Params::TM] = static_cast<RealT>(0.0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> zero_time_model(
            &bus,
            zero_time_constants);
        success *= (zero_time_model.verify() == 0);
        success *= (bus.allocate() == 0);
        success *= (bus.initialize() == 0);
        success *= (zero_time_model.allocate() == 0);
        success *= (zero_time_model.initialize() == 0);
        success *= (zero_time_model.evaluateResidual() == 0);
        success *= allResidualsZero(zero_time_model);

        return success.report(__func__);
      }

      TestOutcome initializesFromPowerFlowAndPublishesSignals()
      {
        TestStatus success = true;

        const ScalarT vr{0.8};
        const ScalarT vi{0.6};
        const ScalarT p0{0.4};
        const ScalarT q0{-0.1};
        const RealT   mva{50.0};

        PhasorDynamics::Bus<ScalarT, IdxT> bus(vr, vi);
        bus.allocate();
        bus.initialize();

        auto data                    = makeData();
        data.parameters[Params::mva] = mva;
        data.parameters[Params::p0]  = static_cast<RealT>(p0);
        data.parameters[Params::q0]  = static_cast<RealT>(q0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        ScalarT ipcmd_value{-1.0};
        ScalarT iqcmd_value{1.0};
        IdxT    ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT    iqcmd_index = ipcmd_index + 1;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pbranch_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qbranch_node;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);

        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);
        regca.getSignals().template assignSignalNode<Vars::IR>(&ir_node);
        regca.getSignals().template assignSignalNode<Vars::II>(&ii_node);
        regca.getSignals().template assignSignalNode<Vars::PBR>(&pbranch_node);
        regca.getSignals().template assignSignalNode<Vars::QBR>(&qbranch_node);

        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);
        success *= (regca.evaluateResidual() == 0);

        // Fixed answer key for Vr = 0.8, Vi = 0.6, P0 = 0.4, Q0 = -0.1,
        // and a 50 MVA component base on the 100 MVA system base.
        const auto* y  = regca.y().getData();
        success       *= scalarMatches(y[index(Vars::VM)], static_cast<ScalarT>(1.0), "VM");
        success       *= scalarMatches(y[index(Vars::VT)], static_cast<ScalarT>(1.0), "VT");
        success       *= scalarMatches(y[index(Vars::IP)], static_cast<ScalarT>(0.8), "IP");
        success       *= scalarMatches(y[index(Vars::IQ)], static_cast<ScalarT>(-0.2), "IQ");
        success       *= scalarMatches(y[index(Vars::IR)], static_cast<ScalarT>(0.26), "IR");
        success       *= scalarMatches(y[index(Vars::II)], static_cast<ScalarT>(0.32), "II");
        success       *= scalarMatches(y[index(Vars::PBR)], p0, "PBR");
        success       *= scalarMatches(y[index(Vars::QBR)], q0, "QBR");
        success       *= scalarMatches(ipcmd_node.read(), p0, "ipcmd signal");
        success       *= scalarMatches(iqcmd_node.read(), q0, "iqcmd signal");
        success       *= scalarMatches(ir_node.read(), static_cast<ScalarT>(0.26), "ir signal");
        success       *= scalarMatches(ii_node.read(), static_cast<ScalarT>(0.32), "ii signal");
        success       *= scalarMatches(pbranch_node.read(), p0, "pbranch signal");
        success       *= scalarMatches(qbranch_node.read(), q0, "qbranch signal");

        success *= allResidualsZero(regca);
        return success.report(__func__);
      }

      TestOutcome unconnectedCommandsRemainConstant()
      {
        TestStatus success = true;

        auto data                   = makeData();
        data.parameters[Params::p0] = static_cast<RealT>(0.6);
        data.parameters[Params::q0] = static_cast<RealT>(0.2);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        success       *= (regca.allocate() == 0);
        success       *= (regca.verify() == 0);
        success       *= (regca.initialize() == 0);
        success       *= (regca.evaluateResidual() == 0);
        const auto* y  = regca.y().getData();

        success *= isEqual(y[index(Vars::IP)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(y[index(Vars::IQ)], static_cast<ScalarT>(0.2), kTol);
        success *= allResidualsZero(regca);

        return success.report(__func__);
      }

      TestOutcome initializesNearHighVoltageLimit()
      {
        TestStatus success = true;

        // Half the transition margin keeps the root distinct from the margin, so
        // a swapped-argument HVRCM residual cannot satisfy this answer key.
        const RealT   voltage_margin = HALF<RealT> * kHvrcmTransition;
        const ScalarT terminal_voltage{kHvrcmVoltageLimit - voltage_margin};
        const ScalarT expected_extra_current{kHvrcmHalfTransitionCurrent};
        const ScalarT q0{0.1};

        auto data                   = makeData();
        data.parameters[Params::q0] = static_cast<RealT>(q0);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(terminal_voltage, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);
        success *= (regca.evaluateResidual() == 0);
        success *= allResidualsZero(regca);

        const auto* y  = regca.y().getData();
        success       *= scalarMatches(y[index(Vars::IQEXTRA)],
                                 expected_extra_current,
                                 "IQEXTRA smooth-constraint root");
        success       *= scalarMatches(y[index(Vars::IQ)],
                                 q0 / terminal_voltage + expected_extra_current,
                                 "IQ preserves Q0 after HVRCM compensation");
        success       *= scalarMatches(y[index(Vars::QBR)], q0, "QBR");

        return success.report(__func__);
      }

      TestOutcome rejectsInitializationAtOrAboveHighVoltageLimit()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing initialization at and above Vhvmax. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        const std::array<RealT, 2> terminal_voltages{{
            kHvrcmVoltageLimit,
            kHvrcmVoltageLimit + static_cast<RealT>(0.1),
        }};

        for (const RealT terminal_voltage : terminal_voltages)
        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(terminal_voltage, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeData());
          success *= (regca.allocate() == 0);
          success *= (regca.verify() == 0);

          auto* y = regca.y().getData();
          std::fill(y,
                    y + static_cast<size_t>(regca.size()),
                    static_cast<ScalarT>(0.25));
          regca.y().setDataUpdated();

          success *= (regca.initialize() != 0);
          for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
          {
            success *= isEqual(y[i], static_cast<ScalarT>(0.25), kTol);
          }
        }

        return success.report(__func__);
      }

      // Verifies initialization rejects an operating point below VA0. VA0 is
      // 0.4 and VA1 is 0.9 in makeData.
      TestOutcome rejectsInitializationBelowLvacmBreakpoint()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing initialization below the LVACM breakpoint. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        auto data                   = makeData();
        data.parameters[Params::q0] = static_cast<RealT>(0.1);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.2, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() > 0);

        return success.report(__func__);
      }

      // Verifies initialization rejects an operating point on the LVACM ramp.
      // VA0 is 0.4 and VA1 is 0.9 in makeData.
      TestOutcome rejectsInitializationWithActiveLvacm()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing initialization with active LVACM. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        const ScalarT p0{0.13};

        auto data                   = makeData();
        data.parameters[Params::p0] = static_cast<RealT>(p0);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.65, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() > 0);

        return success.report(__func__);
      }

      // Verifies VA1 itself is an admissible initialization boundary and the
      // resulting operating point reproduces the P0/Q0 injection.
      TestOutcome initializesAtLvacmUpperBreakpoint()
      {
        TestStatus success = true;

        const ScalarT p0{0.2};
        const ScalarT q0{0.1};

        auto data                   = makeData();
        data.parameters[Params::p0] = static_cast<RealT>(p0);
        data.parameters[Params::q0] = static_cast<RealT>(q0);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.9, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);
        success *= (regca.evaluateResidual() == 0);
        success *= allResidualsZero(regca);

        const auto* y  = regca.y().getData();
        success       *= scalarMatches(y[index(Vars::PBR)], p0, "PBR holds P0");
        success       *= scalarMatches(y[index(Vars::QBR)], q0, "QBR holds Q0");

        return success.report(__func__);
      }

      // Verifies a degenerate terminal voltage is rejected before any state is
      // written, since every current command divides by the voltage magnitude.
      TestOutcome rejectsZeroTerminalVoltage()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing that a zero terminal voltage is rejected. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeData());
        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() > 0);

        return success.report(__func__);
      }

      TestOutcome signalVerification()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing REGCA command signal verification. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.25};
        ScalarT                                   iqcmd_value{-0.10};
        IdxT                                      ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT                                      iqcmd_index = ipcmd_index + 1;

        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        success *= (regca.verify() > 0);

        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        success *= (regca.verify() == 0);

        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);
        success *= (regca.verify() > 0);

        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        success *= (regca.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome nullBusVerification()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing REGCA verification without a bus. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(previous_verbosity);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(nullptr, makeData());
        success *= (regca.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome busInjectionUsesSystemBase()
      {
        TestStatus success = true;

        const ScalarT vr{0.8};
        const ScalarT vi{0.6};
        const ScalarT p0{0.4};
        const ScalarT q0{-0.1};
        const RealT   mva{50.0};

        PhasorDynamics::Bus<ScalarT, IdxT> bus(vr, vi);
        bus.allocate();
        bus.initialize();

        auto data                    = makeData();
        data.parameters[Params::mva] = mva;
        data.parameters[Params::p0]  = static_cast<RealT>(p0);
        data.parameters[Params::q0]  = static_cast<RealT>(q0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.initialize() == 0);

        bus.evaluateResidual();
        success       *= (regca.evaluateResidual() == 0);
        const auto* y  = regca.y().getData();

        success *= isEqual(bus.Ir(), y[index(Vars::IR)], kTol);
        success *= isEqual(bus.Ii(), y[index(Vars::II)], kTol);
        success *= isEqual(bus.Ir(), static_cast<ScalarT>(0.26), kTol);
        success *= isEqual(bus.Ii(), static_cast<ScalarT>(0.32), kTol);
        success *= allResidualsZero(regca);

        return success.report(__func__);
      }

      TestOutcome residualEquations()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(kStateVr, kStateVi);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeDynamicData());

        ScalarT ipcmd_value{0.0};
        ScalarT iqcmd_value{0.0};
        IdxT    ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT    iqcmd_index = ipcmd_index + 1;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);

        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);

        // initialize() seeds the attached command nodes from p0 and q0. Drive
        // them to the values the answer key assumes.
        ipcmd_value = static_cast<ScalarT>(kStateIpcmd);
        iqcmd_value = static_cast<ScalarT>(kStateIqcmd);

        setResidualState(regca);
        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);

        // Answer key for the state set by setResidualState() with the commands
        // above. Each entry is the expected value of the named README equation.
        const std::vector<ScalarT> res_answer = {
            0.865,  // f[0]:  -VM' + (VT - VM) / TM
            1.52,   // f[1]:  -IQ' + max(fq, Rqmin)
            0.22,   // f[2]:  -IP' + clamp(fp, LP, UP)
            -0.035, // f[3]:  -VT^2 + Vr^2 + Vi^2
            0.25,   // f[4]:  -VT*IR + Vi*(IQ - IQEXTRA) + Vr*IP*linseg(VT)
            0.251,  // f[5]:  -VT*II - Vr*(IQ - IQEXTRA) + Vi*IP*linseg(VT)
            -0.03,  // f[6]:  smooth HVRCM constraint
            0.35,   // f[7]:  -IL + linseg(VM, VL0, VL1, IL1)
            -69.6,  // f[8]:  -LP - Rpmax - (Mp - Rpmax)*sigmoid(IP)
            -0.5,   // f[9]:  -UP + Mp*(1 - sigmoid(IP)) + Rpmax*sigmoid(IP)*sigmoid(IL - IP)
            0.0,    // f[10]: -PBR + Vr*IR + Vi*II
            0.0,    // f[11]: -QBR + Vi*IR - Vr*II
        };

        const auto& residual  = regca.getResidual();
        success              *= (static_cast<size_t>(residual.getSize()) == res_answer.size());

        const auto* f = residual.getData();
        for (size_t i = 0; i < res_answer.size(); ++i)
        {
          if (!isEqual(f[i], res_answer[i], kTol))
          {
            std::cout << "Incorrect result for residual " << i << ": "
                      << std::setprecision(15) << f[i] << " != " << res_answer[i]
                      << "\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

      TestOutcome highVoltageReactiveCurrentConstraint()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeData());
        success *= (regca.allocate() == 0);
        success *= (regca.initialize() == 0);

        struct HvrcmCase
        {
          RealT voltage;
          RealT extra_current;
          bool  is_root;
        };

        const std::array<HvrcmCase, 6> cases{{
            {kHvrcmVoltageLimit - static_cast<RealT>(0.1), ZERO<RealT>, true},
            {kHvrcmVoltageLimit - kHvrcmTransition, kHvrcmTransition, true},
            {kHvrcmVoltageLimit - static_cast<RealT>(0.1), static_cast<RealT>(0.05), false},
            {kHvrcmVoltageLimit - static_cast<RealT>(0.1), static_cast<RealT>(-0.05), false},
            {kHvrcmVoltageLimit, ZERO<RealT>, false},
            {kHvrcmVoltageLimit + static_cast<RealT>(0.1), static_cast<RealT>(0.1), false},
        }};

        auto* y = regca.y().getData();
        for (const auto& test_case : cases)
        {
          y[index(Vars::VT)]      = static_cast<ScalarT>(test_case.voltage);
          y[index(Vars::IQEXTRA)] = static_cast<ScalarT>(test_case.extra_current);
          regca.y().setDataUpdated();

          bus.evaluateResidual();
          success *= (regca.evaluateResidual() == 0);

          const auto* f = regca.getResidual().getData();
          const bool  at_root =
              isEqual(f[index(Vars::IQEXTRA)], static_cast<ScalarT>(0.0), kTol);
          success *= (at_root == test_case.is_root);
        }

        return success.report(__func__);
      }

      TestOutcome highVoltageReactiveCurrentJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        TestStatus success = true;

        PhasorDynamics::Bus<DepVar, IdxT> bus(DepVar{1.0}, DepVar{0.0});
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<DepVar, IdxT> regca(&bus, makeData());
        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);

        const size_t voltage_index = index(Vars::VT);
        const size_t current_index = index(Vars::IQEXTRA);
        auto*        y             = regca.y().getData();
        y[voltage_index].setValue(kHvrcmVoltageLimit - kHvrcmTransition);
        y[voltage_index].setVariableNumber(voltage_index);
        y[current_index].setValue(kHvrcmTransition);
        y[current_index].setVariableNumber(current_index);
        regca.y().setDataUpdated();

        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);
        const auto& dependencies =
            regca.getResidual().getData()[current_index].getDependencies();

        success *= isEqual(dependencies.at(voltage_index),
                           static_cast<RealT>(0.5),
                           static_cast<RealT>(1.0e-10));
        success *= isEqual(dependencies.at(current_index),
                           static_cast<RealT>(-0.5),
                           static_cast<RealT>(1.0e-10));

        return success.report(__func__);
      }

      // Positive initial reactive power selects the upper IQ recovery-rate limit.
      TestOutcome positiveInitialReactivePowerSelectsUpperIqRateLimit()
      {
        TestStatus success = true;

        auto data                   = makeDynamicData();
        data.parameters[Params::q0] = static_cast<RealT>(0.2);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        ScalarT                                   iqcmd_value{0.2};
        IdxT                                      iqcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);

        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);

        iqcmd_value = static_cast<ScalarT>(0.4);
        bus.evaluateResidual();
        success       *= (regca.evaluateResidual() == 0);
        const auto* f  = regca.getResidual().getData();

        success *= scalarMatches(f[index(Vars::IQ)],
                                 static_cast<ScalarT>(0.5),
                                 "upper IQ rate limit");

        return success.report(__func__);
      }

      // Disabling LVPL removes IL from the active-current upper rate bound.
      TestOutcome disabledLvplRemovesIlDependence()
      {
        TestStatus success = true;

        auto data                   = makeDynamicData();
        data.parameters[Params::sL] = false;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.initialize() == 0);

        auto* y            = regca.y().getData();
        y[index(Vars::IP)] = static_cast<ScalarT>(0.5);
        y[index(Vars::IL)] = static_cast<ScalarT>(0.1);
        y[index(Vars::UP)] = static_cast<ScalarT>(0.0);
        regca.y().setDataUpdated();

        bus.evaluateResidual();
        success       *= (regca.evaluateResidual() == 0);
        const auto* f  = regca.getResidual().getData();

        success *= scalarMatches(f[index(Vars::UP)],
                                 static_cast<ScalarT>(0.7),
                                 "UP target with LVPL disabled");

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto dependency_tracking_jacobian = dependencyTrackingJacobian();
        auto enzyme_jacobian              = enzymeJacobian();

        success          *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
        const auto nrows  = std::min(dependency_tracking_jacobian.size(), enzyme_jacobian.size());

        for (size_t i = 0; i < nrows; ++i)
        {
          success *= isEqual(dependency_tracking_jacobian[i],
                             enzyme_jacobian[i],
                             static_cast<RealT>(kTol));
        }

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Converter::RegcaParameters;
      using Vars   = PhasorDynamics::Converter::RegcaInternalVariables;
      using Ext    = PhasorDynamics::Converter::RegcaExternalVariables;

      static size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      // Operating point shared by the residual answer key and both Jacobian
      // paths, so the scalar and dependency-tracking states cannot drift apart.
      static constexpr RealT kStateVr    = static_cast<RealT>(0.95);
      static constexpr RealT kStateVi    = static_cast<RealT>(0.25);
      static constexpr RealT kStateIpcmd = static_cast<RealT>(0.9);
      static constexpr RealT kStateIqcmd = static_cast<RealT>(0.1);

      /// Vhvmax in makeData() and in makeDynamicData().
      static constexpr RealT kHvrcmVoltageLimit        = static_cast<RealT>(1.2);
      static constexpr RealT kDynamicHvrcmVoltageLimit = static_cast<RealT>(1.3);

      /// log(2) / Math::MU. At this voltage margin the smooth HVRCM constraint
      /// root equals the margin itself.
      static constexpr RealT kHvrcmTransition = static_cast<RealT>(0.0028881132523331052);

      /// Constraint root at half the transition margin, where root and margin
      /// differ. Held as a literal so the tests stay free of <cmath>.
      static constexpr RealT kHvrcmHalfTransitionCurrent =
          static_cast<RealT>(0.005116446572081316);

      /// HVRCM transition operating point on the makeDynamicData() limit.
      static constexpr RealT kHvrcmTransitionVoltage =
          kDynamicHvrcmVoltageLimit - kHvrcmTransition;

      /// Internal variables in `RegcaInternalVariables` order.
      static constexpr std::array<RealT, static_cast<size_t>(Vars::MAXIMUM)> kStateY = {
          static_cast<RealT>(0.65),   // VM
          static_cast<RealT>(-0.2),   // IQ
          static_cast<RealT>(0.85),   // IP
          static_cast<RealT>(1.0),    // VT
          static_cast<RealT>(0.5),    // IR
          static_cast<RealT>(0.18),   // II
          static_cast<RealT>(0.03),   // IQEXTRA
          static_cast<RealT>(0.2),    // IL
          static_cast<RealT>(-0.4),   // LP
          static_cast<RealT>(0.5),    // UP
          static_cast<RealT>(0.52),   // PBR
          static_cast<RealT>(-0.046), // QBR
      };

      /// Derivatives of the three differential states, which lead the enum.
      static constexpr std::array<RealT, 3> kStateYp = {
          static_cast<RealT>(0.01),  // VM
          static_cast<RealT>(-0.02), // IQ
          static_cast<RealT>(0.03),  // IP
      };

      auto makeData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        using Buses = PhasorDynamics::Converter::RegcaBuses;
        using Mon   = PhasorDynamics::Converter::RegcaMonitorableVariables;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> data;
        data.device_class          = "Regca";
        data.disambiguation_string = "regca_test";
        data.buses[Buses::bus]     = 1;
        data.monitored_variables.insert(Mon::ir);
        data.monitored_variables.insert(Mon::ii);
        data.monitored_variables.insert(Mon::p);
        data.monitored_variables.insert(Mon::q);

        data.parameters[Params::p0]     = static_cast<RealT>(0.0);
        data.parameters[Params::q0]     = static_cast<RealT>(0.0);
        data.parameters[Params::mva]    = static_cast<RealT>(100.0);
        data.parameters[Params::Tg]     = static_cast<RealT>(0.02);
        data.parameters[Params::TM]     = static_cast<RealT>(0.02);
        data.parameters[Params::Rqmax]  = static_cast<RealT>(999.0);
        data.parameters[Params::Rqmin]  = static_cast<RealT>(-999.0);
        data.parameters[Params::Rpmax]  = static_cast<RealT>(999.0);
        data.parameters[Params::sL]     = true;
        data.parameters[Params::IL1]    = static_cast<RealT>(1.1);
        data.parameters[Params::VL0]    = static_cast<RealT>(0.4);
        data.parameters[Params::VL1]    = static_cast<RealT>(0.9);
        data.parameters[Params::VA0]    = static_cast<RealT>(0.4);
        data.parameters[Params::VA1]    = static_cast<RealT>(0.9);
        data.parameters[Params::Vhvmax] = kHvrcmVoltageLimit;

        return data;
      }

      auto makeDynamicData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        auto data                       = makeData();
        data.parameters[Params::p0]     = static_cast<RealT>(0.6);
        data.parameters[Params::q0]     = static_cast<RealT>(-0.2);
        data.parameters[Params::Tg]     = static_cast<RealT>(0.2);
        data.parameters[Params::TM]     = static_cast<RealT>(0.4);
        data.parameters[Params::Rqmax]  = static_cast<RealT>(0.5);
        data.parameters[Params::Rqmin]  = static_cast<RealT>(-0.6);
        data.parameters[Params::Rpmax]  = static_cast<RealT>(0.7);
        data.parameters[Params::Vhvmax] = kDynamicHvrcmVoltageLimit;
        return data;
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
          if (!isEqual(f[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "REGCA residual row " << i << " is " << f[i] << "\n";
            success = false;
          }
          if (!isEqual(yp[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "REGCA derivative row " << i << " is " << yp[i] << "\n";
            success = false;
          }
        }
        return success;
      }

      bool scalarMatches(ScalarT     actual,
                         ScalarT     expected,
                         const char* label,
                         ScalarT     tol = kTol) const
      {
        if (isEqual(actual, expected, tol))
        {
          return true;
        }

        std::cout << label << " mismatch: " << actual << " != " << expected << "\n";
        return false;
      }

      void setResidualState(PhasorDynamics::Converter::Regca<ScalarT, IdxT>& regca)
      {
        auto* y  = regca.y().getData();
        auto* yp = regca.yp().getData();

        for (size_t i = 0; i < kStateY.size(); ++i)
        {
          y[i] = static_cast<ScalarT>(kStateY[i]);
        }
        for (size_t i = 0; i < kStateYp.size(); ++i)
        {
          yp[i] = static_cast<ScalarT>(kStateYp[i]);
        }

        regca.y().setDataUpdated();
        regca.yp().setDataUpdated();
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      void setResidualStateDep(
          PhasorDynamics::Converter::Regca<DependencyTracking::Variable, IdxT>& regca,
          PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>&              bus)
      {
        auto* bus_y = bus.y().getData();
        auto* y     = regca.y().getData();
        auto* yp    = regca.yp().getData();

        bus_y[0].setValue(kStateVr);
        bus_y[1].setValue(kStateVi);

        for (size_t i = 0; i < kStateY.size(); ++i)
        {
          y[i].setValue(kStateY[i]);
        }
        for (size_t i = 0; i < kStateYp.size(); ++i)
        {
          yp[i].setValue(kStateYp[i]);
        }

        bus.y().setDataUpdated();
        regca.y().setDataUpdated();
        regca.yp().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        auto data                    = makeDynamicData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        PhasorDynamics::Bus<DepVar, IdxT>              bus(DepVar{kStateVr}, DepVar{kStateVi});
        PhasorDynamics::Converter::Regca<DepVar, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<DepVar, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> iqcmd_node;
        DepVar                                   ipcmd_value{kStateIpcmd};
        DepVar                                   iqcmd_value{kStateIqcmd};
        IdxT                                     ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT                                     iqcmd_index = ipcmd_index + 1;

        bus.allocate();
        regca.allocate();
        bus.initialize();
        setResidualStateDep(regca, bus);
        auto* y     = regca.y().getData();
        auto* yp    = regca.yp().getData();
        auto* bus_y = bus.y().getData();

        y[index(Vars::VT)].setValue(kHvrcmTransitionVoltage);
        y[index(Vars::IQEXTRA)].setValue(kHvrcmTransition);

        for (IdxT i = 0; i < regca.size(); ++i)
        {
          y[static_cast<size_t>(i)].setVariableNumber(i);
          yp[static_cast<size_t>(i)].setVariableNumber(i);
        }
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus_y[static_cast<size_t>(i)].setVariableNumber(i + regca.size());
        }
        regca.y().setDataUpdated();
        regca.yp().setDataUpdated();
        bus.y().setDataUpdated();
        ipcmd_value.setVariableNumber(ipcmd_index);
        iqcmd_value.setVariableNumber(iqcmd_index);

        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        bus.evaluateResidual();
        regca.evaluateResidual();

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(
            static_cast<size_t>(regca.size() + bus.size()));
        const auto* f = regca.getResidual().getData();
        for (IdxT i = 0; i < regca.size(); ++i)
        {
          dependencies[static_cast<size_t>(i)] =
              f[static_cast<size_t>(i)].getDependencies();
        }
        dependencies[static_cast<size_t>(regca.size())]     = bus.Ir().getDependencies();
        dependencies[static_cast<size_t>(regca.size() + 1)] = bus.Ii().getDependencies();

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian()
      {
        auto data                    = makeDynamicData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(kStateVr, kStateVi);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{kStateIpcmd};
        ScalarT                                   iqcmd_value{kStateIqcmd};
        IdxT                                      ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT                                      iqcmd_index = ipcmd_index + 1;

        bus.allocate();
        regca.allocate();
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + regca.size());
          bus.setResidualIndex(i, i + regca.size());
        }

        bus.initialize();
        setResidualState(regca);
        auto* y                 = regca.y().getData();
        y[index(Vars::VT)]      = static_cast<ScalarT>(kHvrcmTransitionVoltage);
        y[index(Vars::IQEXTRA)] = static_cast<ScalarT>(kHvrcmTransition);
        regca.y().setDataUpdated();
        regca.updateTime(0.0, 1.0);

        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        bus.evaluateResidual();
        regca.evaluateResidual();
        regca.evaluateJacobian();

        regca.constructCsr();
        return MapFromCsr(regca.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
