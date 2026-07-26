#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
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

      // Tolerance for values the model computes without a smooth approximation.
      static constexpr ScalarT kTol = std::numeric_limits<ScalarT>::epsilon();

      // initialize() divides by linseg(VT, VA0, VA1, 1). Above VA1 that value is
      // short of one by 3e-13, and the residue reaches every state seeded from
      // P0. Comparisons against exact power-flow values carry it.
      static constexpr ScalarT kInitTol = static_cast<ScalarT>(1.0e-9);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> minimal(&bus);
        success *= (minimal.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (minimal.getMonitor() == nullptr);
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

        auto missing = makeData();
        missing.parameters.erase(Params::Tg);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing_model(&bus, missing);
        success *= (missing_model.verify() > 0);

        auto missing_p0 = makeData();
        missing_p0.parameters.erase(Params::P0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing_p0_model(&bus, missing_p0);
        success *= (missing_p0_model.verify() > 0);

        auto missing_q0 = makeData();
        missing_q0.parameters.erase(Params::Q0);
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
        data.parameters[Params::P0]  = static_cast<RealT>(p0);
        data.parameters[Params::Q0]  = static_cast<RealT>(q0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        ScalarT ipcmd_value{-1.0};
        ScalarT iqcmd_value{1.0};
        IdxT    ipcmd_index = 20;
        IdxT    iqcmd_index = 21;

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

        const ScalarT vt = std::sqrt(vr * vr + vi * vi);
        const ScalarT lvacm =
            Math::linseg(vt, static_cast<RealT>(0.4), static_cast<RealT>(0.9), ONE<RealT>);
        const ScalarT ipcmd0       = toComponentBase(p0 / vt, mva) / lvacm;
        const ScalarT iqcmd0       = toComponentBase(q0 / vt, mva);
        const ScalarT qnet0        = iqcmd0;
        const ScalarT ir0          = (vr * ipcmd0 * lvacm + vi * qnet0) / vt;
        const ScalarT ii0          = (vi * ipcmd0 * lvacm - vr * qnet0) / vt;
        const ScalarT ipcmd_signal = toSystemBase(ipcmd0, mva);
        const ScalarT iqcmd_signal = toSystemBase(iqcmd0, mva);
        const ScalarT ir_signal    = toSystemBase(ir0, mva);
        const ScalarT ii_signal    = toSystemBase(ii0, mva);
        const auto*   y            = regca.y().getData();

        success *= scalarMatches(y[index(Vars::VM)], vt, "VM");
        success *= scalarMatches(y[index(Vars::VT)], vt, "VT");
        success *= scalarMatches(y[index(Vars::IP)], ipcmd0, "IP");
        success *= scalarMatches(y[index(Vars::IQ)], iqcmd0, "IQ");
        success *= scalarMatches(y[index(Vars::IR)], ir_signal, "IR");
        success *= scalarMatches(y[index(Vars::II)], ii_signal, "II");
        success *= scalarMatches(y[index(Vars::PBR)], p0, "PBR");
        success *= scalarMatches(y[index(Vars::QBR)], q0, "QBR");
        success *= scalarMatches(ipcmd_node.read(), ipcmd_signal, "ipcmd signal");
        success *= scalarMatches(iqcmd_node.read(), iqcmd_signal, "iqcmd signal");
        success *= scalarMatches(ir_node.read(), ir_signal, "ir signal");
        success *= scalarMatches(ii_node.read(), ii_signal, "ii signal");
        success *= scalarMatches(pbranch_node.read(), p0, "pbranch signal");
        success *= scalarMatches(qbranch_node.read(), q0, "qbranch signal");

        success *= allResidualsZero(regca);
        return success.report(__func__);
      }

      TestOutcome baseSignals()
      {
        TestStatus success = true;

        const ScalarT p0{0.4};
        const ScalarT q0{-0.1};
        const RealT   mva{50.0};

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeData();
        data.parameters[Params::mva] = mva;
        data.parameters[Params::P0]  = static_cast<RealT>(p0);
        data.parameters[Params::Q0]  = static_cast<RealT>(q0);

        ScalarT ipcmd_value{99.0};
        ScalarT iqcmd_value{99.0};
        IdxT    ipcmd_index = 30;
        IdxT    iqcmd_index = 31;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pbranch_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qbranch_node;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);
        regca.getSignals().template assignSignalNode<Vars::IR>(&ir_node);
        regca.getSignals().template assignSignalNode<Vars::II>(&ii_node);
        regca.getSignals().template assignSignalNode<Vars::PBR>(&pbranch_node);
        regca.getSignals().template assignSignalNode<Vars::QBR>(&qbranch_node);

        success       *= (regca.allocate() == 0);
        success       *= (regca.verify() == 0);
        success       *= (regca.initialize() == 0);
        success       *= (regca.evaluateResidual() == 0);
        const auto* y  = regca.y().getData();

        success *= scalarMatches(y[index(Vars::IP)], static_cast<ScalarT>(0.8), "IP", kInitTol);
        success *= scalarMatches(y[index(Vars::IQ)], static_cast<ScalarT>(-0.2), "IQ");
        success *= scalarMatches(y[index(Vars::IR)], static_cast<ScalarT>(0.4), "IR");
        success *= scalarMatches(y[index(Vars::II)], static_cast<ScalarT>(0.1), "II");
        success *= scalarMatches(y[index(Vars::PBR)], p0, "PBR");
        success *= scalarMatches(y[index(Vars::QBR)], q0, "QBR");
        success *= scalarMatches(ipcmd_node.read(), p0, "ipcmd signal", kInitTol);
        success *= scalarMatches(iqcmd_node.read(), q0, "iqcmd signal");
        success *= scalarMatches(ir_node.read(), static_cast<ScalarT>(0.4), "ir signal");
        success *= scalarMatches(ii_node.read(), static_cast<ScalarT>(0.1), "ii signal");
        success *= scalarMatches(pbranch_node.read(), p0, "pbranch signal");
        success *= scalarMatches(qbranch_node.read(), q0, "qbranch signal");
        success *= allResidualsZero(regca);

        return success.report(__func__);
      }

      TestOutcome unconnectedCommandsRemainConstant()
      {
        TestStatus success = true;

        auto data                   = makeData();
        data.parameters[Params::P0] = static_cast<RealT>(0.6);
        data.parameters[Params::Q0] = static_cast<RealT>(0.2);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        success       *= (regca.allocate() == 0);
        success       *= (regca.verify() == 0);
        success       *= (regca.initialize() == 0);
        success       *= (regca.evaluateResidual() == 0);
        const auto* y  = regca.y().getData();

        success *= isEqual(y[index(Vars::IP)], static_cast<ScalarT>(0.6), kInitTol);
        success *= isEqual(y[index(Vars::IQ)], static_cast<ScalarT>(0.2), kTol);
        success *= allResidualsZero(regca);

        return success.report(__func__);
      }

      TestOutcome externalCommandsDriveRuntimeResidual()
      {
        TestStatus success = true;

        auto data                   = makeData();
        data.parameters[Params::P0] = static_cast<RealT>(0.5);
        data.parameters[Params::Q0] = static_cast<RealT>(-0.1);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ScalarT ipcmd_value{0.0};
        ScalarT iqcmd_value{0.0};
        IdxT    ipcmd_index = 22;
        IdxT    iqcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);

        ipcmd_value = static_cast<ScalarT>(0.55);
        iqcmd_value = static_cast<ScalarT>(-0.05);
        bus.evaluateResidual();
        success       *= (regca.evaluateResidual() == 0);
        const auto* f  = regca.getResidual().getData();

        // Both commands move 0.05 per unit and sit far inside the rate limits,
        // so each residual is the unlimited command error over Tg.
        success *= isEqual(f[index(Vars::IP)], static_cast<ScalarT>(2.5), kInitTol);
        success *= isEqual(f[index(Vars::IQ)], static_cast<ScalarT>(2.5), kTol);
        success *= (f[index(Vars::IP)] > ZERO<RealT>);
        success *= (f[index(Vars::IQ)] > ZERO<RealT>);

        return success.report(__func__);
      }

      // Verifies the initial point is residual-consistent above Vhvmax and still
      // delivers the P0/Q0 injection. Vhvmax is 1.2 in makeData.
      TestOutcome initializesAboveHighVoltageLimit()
      {
        TestStatus success = true;

        const ScalarT q0{0.1};

        auto data                   = makeData();
        data.parameters[Params::Q0] = static_cast<RealT>(q0);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.3, 0.0);
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
                                 Math::ramp(y[index(Vars::VT)] - static_cast<RealT>(1.2)),
                                 "IQEXTRA seeded from the HVRCM ramp");

        // The command absorbs the HVRCM current, so the injection still matches
        // the power-flow values.
        success *= scalarMatches(y[index(Vars::PBR)], static_cast<ScalarT>(0.0), "PBR holds P0");
        success *= scalarMatches(y[index(Vars::QBR)], q0, "QBR holds Q0");

        return success.report(__func__);
      }

      // Verifies initialization rejects an operating point below VA0. VA0 is
      // 0.4 and VA1 is 0.9 in makeData.
      TestOutcome rejectsInitializationBelowLvacmBreakpoint()
      {
        TestStatus success = true;

        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing initialization below the LVACM breakpoint. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);

        auto data                   = makeData();
        data.parameters[Params::Q0] = static_cast<RealT>(0.1);

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

        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing initialization with active LVACM. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);

        const ScalarT p0{0.13};

        auto data                   = makeData();
        data.parameters[Params::P0] = static_cast<RealT>(p0);

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
        data.parameters[Params::P0] = static_cast<RealT>(p0);
        data.parameters[Params::Q0] = static_cast<RealT>(q0);

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
        success       *= scalarMatches(y[index(Vars::PBR)], p0, "PBR holds P0", kInitTol);
        success       *= scalarMatches(y[index(Vars::QBR)], q0, "QBR holds Q0", kInitTol);

        return success.report(__func__);
      }

      // Verifies a degenerate terminal voltage is rejected before any state is
      // written, since every current command divides by the voltage magnitude.
      TestOutcome rejectsZeroTerminalVoltage()
      {
        TestStatus success = true;

        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing that a zero terminal voltage is rejected. "
                    << "Logged errors are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);

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

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.25};
        ScalarT                                   iqcmd_value{-0.10};
        IdxT                                      ipcmd_index = 24;
        IdxT                                      iqcmd_index = 25;

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
        data.parameters[Params::P0]  = static_cast<RealT>(p0);
        data.parameters[Params::Q0]  = static_cast<RealT>(q0);

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

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.95, 0.25);
        bus.allocate();
        bus.initialize();

        ScalarT ipcmd_value{0.9};
        ScalarT iqcmd_value{0.1};
        IdxT    ipcmd_index = 26;
        IdxT    iqcmd_index = 27;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeDynamicData());
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        success *= (regca.allocate() == 0);
        success *= (regca.initialize() == 0);

        // initialize() seeds the attached command nodes from P0 and Q0. Drive
        // them to the values the answer key assumes.
        ipcmd_value = static_cast<ScalarT>(0.9);
        iqcmd_value = static_cast<ScalarT>(0.1);

        setResidualState(regca);
        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);

        // Answer key for the state set by setResidualState() with the commands
        // above. Each entry is the exact piecewise value of the named equation
        // in README.md.
        const std::vector<ScalarT> res_answer = {
            0.29,   // f[0]:  -VM' + (VT - VM) / TM
            1.52,   // f[1]:  -IQ' + max(fq, Rqmin)
            0.22,   // f[2]:  -IP' + clamp(fp, LP, UP)
            0.0046, // f[3]:  -VT^2 + Vr^2 + Vi^2
            0.26,   // f[4]:  -VT*IR + Vi*(IQ - IQEXTRA) + Vr*IP*linseg(VT)
            0.2546, // f[5]:  -VT*II - Vr*(IQ - IQEXTRA) + Vi*IP*linseg(VT)
            -0.03,  // f[6]:  -IQEXTRA + ramp(VT - Vhvmax)
            0.292,  // f[7]:  -IL + linseg(VM, VL0, VL1, IL1)
            -69.6,  // f[8]:  -LP - Rpmax - (Mp - Rpmax)*sigmoid(IP)
            -0.3,   // f[9]:  -UP + Mp*(1 - sigmoid(IP)) + Rpmax*sigmoid(IP)*sigmoid(IL - IP)
            0.0,    // f[10]: -PBR + Vr*IR + Vi*II
            0.0,    // f[11]: -QBR + Vi*IR - Vr*II
        };

        // Rows 2, 4, 5, 7, and 9 evaluate a CommonMath approximation within one
        // smoothing width of a breakpoint, so they miss the exact piecewise
        // answer. Row 7 is the largest at 4.8e-7.
        const auto                 smooth_tol = static_cast<ScalarT>(1.0e-6);
        const std::vector<ScalarT> res_tol    = {kTol,
                                                 kTol,
                                                 smooth_tol,
                                                 kTol,
                                                 smooth_tol,
                                                 smooth_tol,
                                                 kTol,
                                                 smooth_tol,
                                                 kTol,
                                                 smooth_tol,
                                                 kTol,
                                                 kTol};

        const auto& residual  = regca.getResidual();
        success              *= (static_cast<size_t>(residual.getSize()) == res_answer.size());

        const auto* f = residual.getData();
        for (size_t i = 0; i < res_answer.size(); ++i)
        {
          if (!isEqual(f[i], res_answer[i], res_tol[i]))
          {
            std::cout << "Incorrect result for residual " << i << ": "
                      << std::setprecision(15) << f[i] << " != " << res_answer[i]
                      << "\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

      TestOutcome highVoltageReactiveCurrentRoot()
      {
        TestStatus success = true;

        auto data = makeDynamicData();

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.initialize() == 0);

        auto* y                 = regca.y().getData();
        y[index(Vars::VT)]      = static_cast<ScalarT>(1.35);
        y[index(Vars::IQEXTRA)] = Math::ramp(y[index(Vars::VT)]
                                             - static_cast<RealT>(1.3));
        regca.y().setDataUpdated();

        bus.evaluateResidual();
        success       *= (regca.evaluateResidual() == 0);
        const auto* f  = regca.getResidual().getData();
        success       *= scalarMatches(f[index(Vars::IQEXTRA)],
                                 static_cast<ScalarT>(0.0),
                                 "HVRCM residual root");

        return success.report(__func__);
      }

      TestOutcome limiterBranchCoverage()
      {
        TestStatus success = true;

        {
          auto data                   = makeDynamicData();
          data.parameters[Params::Q0] = static_cast<RealT>(0.2);

          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
          bus.allocate();
          bus.initialize();

          ScalarT                                   iqcmd_value{0.2};
          IdxT                                      iqcmd_index = 28;
          PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
          iqcmd_node.set(&iqcmd_value, &iqcmd_index);

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
          regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

          success *= (regca.allocate() == 0);
          success *= (regca.verify() == 0);
          success *= (regca.initialize() == 0);

          iqcmd_value = static_cast<ScalarT>(0.4);
          bus.evaluateResidual();
          success       *= (regca.evaluateResidual() == 0);
          const auto* y  = regca.y().getData();
          const auto* f  = regca.getResidual().getData();

          const ScalarT tg       = static_cast<RealT>(0.2);
          const ScalarT fq       = (iqcmd_value - y[index(Vars::IQ)]) / tg;
          const ScalarT expected = Math::min(fq, static_cast<RealT>(0.5));

          success *= scalarMatches(f[index(Vars::IQ)],
                                   expected,
                                   "positive IQ branch");
        }

        {
          auto data                   = makeDynamicData();
          data.parameters[Params::sL] = false;

          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
          success *= (regca.allocate() == 0);
          success *= (regca.initialize() == 0);

          bus.evaluateResidual();
          success       *= (regca.evaluateResidual() == 0);
          const auto* y  = regca.y().getData();
          const auto* f  = regca.getResidual().getData();

          const ScalarT ip       = y[index(Vars::IP)];
          const ScalarT sigma_ip = Math::sigmoid(ip);
          const RealT   rpmax    = static_cast<RealT>(0.7);
          const RealT   mp       = static_cast<RealT>(100.0) * rpmax;
          const ScalarT expected = mp * (ONE<RealT> - sigma_ip) + rpmax * sigma_ip;

          success *= scalarMatches(y[index(Vars::UP)], expected, "UP without LVPL");
          success *= scalarMatches(f[index(Vars::UP)],
                                   static_cast<ScalarT>(0.0),
                                   "UP residual without LVPL");
        }

        return success.report(__func__);
      }

      TestOutcome jsonParseAndSystemAssembly()
      {
        TestStatus success = true;

        std::istringstream input(R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "REGCA full model",
    "case_description": "REGCA parser behavior test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    {
      "number": 1,
      "class": "bus",
      "name": "Bus 1",
      "init": { "Vr": 1.0, "Vi": 0.0 },
      "params": { "kv": 1.0 }
    }
  ],
  "devices": [
    {
      "class": "Regca",
      "ports": { "bus": 1 },
      "id": "CV1",
      "params": {
        "P0": 0.0,
        "Q0": 0.0,
        "mva": 100,
        "Tg": 0.02,
        "TM": 0.02,
        "Rqmax": 999.0,
        "Rqmin": -999.0,
        "Rpmax": 999.0,
        "sL": true,
        "IL1": 1.1,
        "VL0": 0.4,
        "VL1": 0.9,
        "VA0": 0.4,
        "VA1": 0.9,
        "Vhvmax": 1.2
      },
      "mon": ["ir", "ii", "p", "q"]
    }
  ]
}
)json");

        auto data               = PhasorDynamics::parseSystemModelData(input);
        success                *= (data.regca.size() == 1);
        const auto& regca_data  = data.regca[0];
        success                *= (regca_data.device_class == "Regca");
        success                *= (regca_data.buses.at(PhasorDynamics::Converter::RegcaBuses::bus)
                    == 1);
        success                *= regca_data.signal_inputs.empty();
        success                *= regca_data.signal_outputs.empty();
        success                *= (std::get_if<double>(&regca_data.parameters.at(Params::P0))
                    != nullptr);
        success                *= (std::get_if<double>(&regca_data.parameters.at(Params::Q0))
                    != nullptr);
        success                *= (std::get_if<size_t>(&regca_data.parameters.at(Params::mva))
                    != nullptr);
        success                *= (std::get_if<bool>(&regca_data.parameters.at(Params::sL))
                    != nullptr);

        PhasorDynamics::SystemModel<double, size_t> system(data);
        success              *= (system.allocate() == 0);
        success              *= (system.initialize() == 0);
        success              *= (system.tagDifferentiable() == 0);
        success              *= (system.evaluateResidual() == 0);
        success              *= (system.evaluateJacobian() == 0);
        success              *= (system.size() == 14);
        const auto* residual  = system.getResidual().getData();
        success              *= isEqual(residual[0], 0.0, static_cast<double>(kTol));
        success              *= isEqual(residual[1], 0.0, static_cast<double>(kTol));
        for (size_t i = 2; i < static_cast<size_t>(system.getResidual().getSize()); ++i)
        {
          success *= isEqual(residual[i], 0.0, static_cast<double>(kTol));
        }

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
                             static_cast<RealT>(1.0e-8));
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

        data.parameters[Params::P0]     = static_cast<RealT>(0.0);
        data.parameters[Params::Q0]     = static_cast<RealT>(0.0);
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
        data.parameters[Params::Vhvmax] = static_cast<RealT>(1.2);

        return data;
      }

      auto makeDynamicData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        auto data                       = makeData();
        data.parameters[Params::P0]     = static_cast<RealT>(0.6);
        data.parameters[Params::Q0]     = static_cast<RealT>(-0.2);
        data.parameters[Params::Tg]     = static_cast<RealT>(0.2);
        data.parameters[Params::TM]     = static_cast<RealT>(0.4);
        data.parameters[Params::Rqmax]  = static_cast<RealT>(0.5);
        data.parameters[Params::Rqmin]  = static_cast<RealT>(-0.6);
        data.parameters[Params::Rpmax]  = static_cast<RealT>(0.7);
        data.parameters[Params::IL1]    = static_cast<RealT>(1.1);
        data.parameters[Params::Vhvmax] = static_cast<RealT>(1.3);
        return data;
      }

      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus, Params param, RealT value)
      {
        auto data              = makeData();
        data.parameters[param] = value;
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> model(&bus, data);
        return model.verify() > 0;
      }

      ScalarT toComponentBase(ScalarT value, RealT mva) const
      {
        return value * static_cast<RealT>(100.0) / mva;
      }

      ScalarT toSystemBase(ScalarT value, RealT mva) const
      {
        return value * mva / static_cast<RealT>(100.0);
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

        y[index(Vars::VM)]      = static_cast<ScalarT>(0.86);
        y[index(Vars::IQ)]      = static_cast<ScalarT>(-0.2);
        y[index(Vars::IP)]      = static_cast<ScalarT>(0.85);
        y[index(Vars::VT)]      = static_cast<ScalarT>(0.98);
        y[index(Vars::IR)]      = static_cast<ScalarT>(0.5);
        y[index(Vars::II)]      = static_cast<ScalarT>(0.18);
        y[index(Vars::IQEXTRA)] = static_cast<ScalarT>(0.03);
        y[index(Vars::IL)]      = static_cast<ScalarT>(0.72);
        y[index(Vars::LP)]      = static_cast<ScalarT>(-0.4);
        y[index(Vars::UP)]      = static_cast<ScalarT>(0.3);
        y[index(Vars::PBR)]     = static_cast<ScalarT>(0.52);
        y[index(Vars::QBR)]     = static_cast<ScalarT>(-0.046);

        yp[index(Vars::VM)] = static_cast<ScalarT>(0.01);
        yp[index(Vars::IQ)] = static_cast<ScalarT>(-0.02);
        yp[index(Vars::IP)] = static_cast<ScalarT>(0.03);

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

        bus_y[0].setValue(0.95);
        bus_y[1].setValue(0.25);

        y[index(Vars::VM)].setValue(0.86);
        y[index(Vars::IQ)].setValue(-0.2);
        y[index(Vars::IP)].setValue(0.85);
        y[index(Vars::VT)].setValue(0.98);
        y[index(Vars::IR)].setValue(0.5);
        y[index(Vars::II)].setValue(0.18);
        y[index(Vars::IQEXTRA)].setValue(0.03);
        y[index(Vars::IL)].setValue(0.72);
        y[index(Vars::LP)].setValue(-0.4);
        y[index(Vars::UP)].setValue(0.3);
        y[index(Vars::PBR)].setValue(0.52);
        y[index(Vars::QBR)].setValue(-0.046);

        yp[index(Vars::VM)].setValue(0.01);
        yp[index(Vars::IQ)].setValue(-0.02);
        yp[index(Vars::IP)].setValue(0.03);

        bus.y().setDataUpdated();
        regca.y().setDataUpdated();
        regca.yp().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        auto data                    = makeDynamicData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        PhasorDynamics::Bus<DepVar, IdxT>              bus(DepVar{0.95}, DepVar{0.25});
        PhasorDynamics::Converter::Regca<DepVar, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<DepVar, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> iqcmd_node;
        DepVar                                   ipcmd_value{0.9};
        DepVar                                   iqcmd_value{0.1};
        IdxT                                     ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT                                     iqcmd_index = static_cast<IdxT>(regca.size() + bus.size() + 1);

        bus.allocate();
        regca.allocate();
        bus.initialize();
        setResidualStateDep(regca, bus);
        auto* y     = regca.y().getData();
        auto* yp    = regca.yp().getData();
        auto* bus_y = bus.y().getData();

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

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(0.95, 0.25);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.9};
        ScalarT                                   iqcmd_value{0.1};
        IdxT                                      ipcmd_index = static_cast<IdxT>(regca.size() + bus.size());
        IdxT                                      iqcmd_index = static_cast<IdxT>(regca.size() + bus.size() + 1);

        bus.allocate();
        regca.allocate();
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + regca.size());
          bus.setResidualIndex(i, i + regca.size());
        }

        bus.initialize();
        setResidualState(regca);
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
