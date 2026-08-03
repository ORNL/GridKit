#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
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

      // Covers accumulated roundoff from smooth inverses, base round trips, and AD evaluation.
      static constexpr RealT kTol = static_cast<RealT>(100) * std::numeric_limits<RealT>::epsilon();

      // A smooth clamp or deadband evaluated exactly at a limit sits log(2)/mu
      // inside it; boundary expectations below are stated with this offset.
      static RealT clampEdgeOffset()
      {
        return std::log(static_cast<RealT>(2)) / Math::MU<RealT>;
      }

      /// Construction, row layout, differential tags, parameters, buses, and signal-link validation use direct contract checks.
      TestOutcome validation()
      {
        TestStatus success = true;
        noteExpectedLogs("Testing invalid REECB configurations. Logged errors and time-constant warnings are expected.");

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(10));
        success *= (empty.getMonitor() == nullptr);

        const std::array<Vars, 10> row_order{{
            Vars::VMEAS,
            Vars::PMEAS,
            Vars::XPIQ,
            Vars::XPIV,
            Vars::QV,
            Vars::PORD,
            Vars::VT,
            Vars::ILMAX,
            Vars::IQCMD,
            Vars::IPCMD,
        }};
        for (size_t row = 0; row < row_order.size(); ++row)
        {
          success *= (index(row_order[row]) == row);
        }

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> configured(&bus, makeData());
        configured.setSystemBase(60.0, 100.0e6);
        success *= (configured.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (configured.getMonitor() != nullptr);
        success *= (configured.verify() == 0);
        success *= (configured.initialize() != 0);
        success *= (configured.allocate() == 0);
        success *= (configured.tagDifferentiable() == 0);

        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          const bool expected = row <= index(Vars::PORD);
          if (configured.tag()[row] != expected)
          {
            std::cout << "REECB differential tag " << row << " mismatch\n";
            success = false;
          }
        }

        const RealT nan = std::numeric_limits<RealT>::quiet_NaN();
        for (const Params parameter : {Params::mva, Params::Trv, Params::Tp, Params::Vref0, Params::Vdip, Params::Vup, Params::dbd1, Params::dbd2, Params::kqv, Params::Iql1, Params::Iqh1, Params::Qmax, Params::Qmin, Params::Kqp, Params::Kqi, Params::Vmax, Params::Vmin, Params::Kvp, Params::Kvi, Params::Tiq, Params::Tpord, Params::dPmax, Params::dPmin, Params::Pmax, Params::Pmin, Params::Imax})
        {
          success *= invalidParameterCase(bus, parameter, nan);
        }

        success *= invalidParameterCase(bus, Params::mva, 0.0);
        success *= invalidParameterCase(bus, Params::Trv, -0.1);
        success *= invalidParameterCase(bus, Params::Vdip, 1.2);
        success *= invalidParameterCase(bus, Params::dbd1, 0.1);
        success *= invalidParameterCase(bus, Params::Iql1, 2.0);
        success *= invalidParameterCase(bus, Params::Qmin, 2.0);
        success *= invalidParameterCase(bus, Params::Vmin, 2.0);
        success *= invalidParameterCase(bus, Params::dPmin, 0.0);
        success *= invalidParameterCase(bus, Params::dPmax, 0.0);
        success *= invalidParameterCase(bus, Params::Pmin, 2.0);
        success *= invalidParameterCase(bus, Params::Imax, 0.0);
        success *= invalidParameterCase(bus, Params::Trv, std::numeric_limits<RealT>::infinity());
        success *= invalidParameterCase(bus, Params::Imax, -std::numeric_limits<RealT>::infinity());
        success *= invalidParameterCase(bus, Params::mva, true);

        // Selectors accept only bool and integer 0/1 encodings.
        for (const Params flag : {Params::PfFlag, Params::VFlag, Params::QFlag, Params::Pqflag})
        {
          success *= !invalidParameterCase(bus, flag, static_cast<IdxT>(0));
          success *= !invalidParameterCase(bus, flag, static_cast<IdxT>(1));
          success *= invalidParameterCase(bus, flag, static_cast<RealT>(0.0));
          success *= invalidParameterCase(bus, flag, static_cast<RealT>(1.0));
          success *= invalidParameterCase(bus, flag, static_cast<IdxT>(2));
          success *= invalidParameterCase(bus, flag, static_cast<RealT>(0.5));
        }

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> busless(nullptr, makeData());
        busless.setSystemBase(60.0, 100.0e6);
        success *= (busless.verify() > 0);
        success *= unlinkedSignalRejected<Ext::PE>(bus);
        success *= unlinkedSignalRejected<Ext::QGEN>(bus);
        success *= unlinkedSignalRejected<Ext::QEXT>(bus);
        success *= unlinkedSignalRejected<Ext::PFAREF>(bus);
        success *= unlinkedSignalRejected<Ext::PREF>(bus);

        auto zero_time                      = makeData();
        zero_time.parameters[Params::Trv]   = 0.0;
        zero_time.parameters[Params::Tp]    = 0.0;
        zero_time.parameters[Params::Tiq]   = 0.0;
        zero_time.parameters[Params::Tpord] = 0.0;
        Fixture<ScalarT> floored(zero_time);
        success *= floored.initialize(0.0, 0.2);
        success *= (floored.evaluate() == 0);
        success *= allResidualsZero(floored.reecb);

        auto* floored_y               = floored.reecb.y().getData();
        floored_y[index(Vars::VMEAS)] = 0.999;
        floored_y[index(Vars::PMEAS)] = 0.199;
        floored_y[index(Vars::QV)]    = 0.001;
        floored_y[index(Vars::PORD)]  = 0.1995;
        floored.reecb.y().setDataUpdated();
        success *= (floored.evaluate() == 0);
        success *= scalarMatches(floored.reecb.getResidual().getData()[index(Vars::VMEAS)], 1.0, "Trv 1 ms floor");
        success *= scalarMatches(floored.reecb.getResidual().getData()[index(Vars::PMEAS)], 1.0, "Tp 1 ms floor");
        success *= scalarMatches(floored.reecb.getResidual().getData()[index(Vars::QV)], -1.0, "Tiq 1 ms floor");
        success *= scalarMatches(floored.reecb.getResidual().getData()[index(Vars::PORD)], 0.5, "Tpord 1 ms floor");

        Data             default_data;
        Fixture<ScalarT> defaulted(default_data);
        success *= defaulted.initialize(0.1, 0.2);
        success *= (defaulted.evaluate() == 0);
        success *= allResidualsZero(defaulted.reecb);
        success *= scalarMatches(defaulted.reecb.y().getData()[index(Vars::ILMAX)], std::sqrt(1.68), "omitted parameter defaults");

        return success.report(__func__);
      }

      /// Nonidentity-base initialization checks known-input preservation, unknown publication, output aliases, latches, and monitor values.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        auto data                    = makeData();
        data.parameters[Params::mva] = 50.0;

        Fixture<ScalarT> fixture(data, 0.8, 0.6);
        fixture.attachAllInputs(99.0);
        fixture.input(Ext::PE)    = 0.3;
        fixture.input(Ext::QGEN)  = -0.05;
        success                  *= fixture.initialize(0.05, 0.25);
        success                  *= (fixture.evaluate() == 0);

        const auto* y  = fixture.reecb.y().getData();
        success       *= scalarMatches(y[index(Vars::VMEAS)], 1.0, "VMEAS");
        success       *= scalarMatches(y[index(Vars::PMEAS)], 0.6, "PMEAS");
        success       *= scalarMatches(y[index(Vars::PORD)], 0.5, "PORD");
        success       *= scalarMatches(y[index(Vars::VT)], 1.0, "VT");
        success       *= scalarMatches(y[index(Vars::ILMAX)], 1.9364916731037085, "ILMAX");
        success       *= (fixture.iqcmd() == 0.05);
        success       *= (fixture.ipcmd() == 0.25);
        success       *= (fixture.input(Ext::PE) == 0.3);
        success       *= (fixture.input(Ext::QGEN) == -0.05);
        success       *= scalarMatches(fixture.input(Ext::QEXT), 0.05, "published QEXT");
        success       *= scalarMatches(fixture.input(Ext::PFAREF), 0.0, "published PFAREF");
        success       *= scalarMatches(fixture.input(Ext::PREF), 0.25, "published PREF");
        success       *= allResidualsZero(fixture.reecb);

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(fixture.reecb.getMonitor());
        std::stringstream output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string header;
        std::string values;
        std::getline(output, header);
        std::getline(output, values);
        success              *= (header == "t,Reecb_reecb_test_iqcmd,Reecb_reecb_test_ipcmd,Reecb_reecb_test_vmeas,Reecb_reecb_test_pmeas");
        const auto monitored  = Tokenizer<RealT>(values, ',')();
        if (monitored.size() == 5)
        {
          success *= scalarMatches(monitored[1], 0.05, "monitored IQCMD");
          success *= scalarMatches(monitored[2], 0.25, "monitored IPCMD");
          success *= scalarMatches(monitored[3], 1.0, "monitored VMEAS");
          success *= scalarMatches(monitored[4], 0.6, "monitored PMEAS");
        }
        else
        {
          std::cout << "REECB monitor emitted " << monitored.size() << " values instead of 5\n";
          success = false;
        }

        Fixture<ScalarT> latched(data, 1.0, 0.0, 100.0e6, false);
        success *= latched.initialize(0.05, 0.25);
        success *= (latched.evaluate() == 0);
        success *= allResidualsZero(latched.reecb);

        auto system_base_data = makeData();
        system_base_data.parameters.erase(Params::mva);
        Fixture<ScalarT> system_base(system_base_data, 1.0, 0.0, 80.0e6);
        success *= system_base.initialize(0.05, 0.25);
        success *= (system_base.evaluate() == 0);
        success *= allResidualsZero(system_base.reecb);
        success *= scalarMatches(system_base.reecb.y().getData()[index(Vars::PMEAS)], 0.25, "omitted mva PMEAS base");
        success *= scalarMatches(system_base.reecb.y().getData()[index(Vars::PORD)], 0.25, "omitted mva PORD base");
        success *= scalarMatches(system_base.reecb.y().getData()[index(Vars::ILMAX)], std::sqrt(3.9375), "omitted mva ILMAX base");

        return success.report(__func__);
      }

      /// Strict inverse-limiter domains, collapsed limits, zero gains, and failure atomicity are checked at accepted and rejected points.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing inadmissible REECB initialization points. Logged errors are expected.");

        success *= initializationRejectedAtomically(makeData(), 0.0, 0.0, 1.0, "zero active-current endpoint");
        success *= initializationRejectedAtomically(makeData(), 0.0, 2.0, 1.0, "zero ILMAX at active-current endpoint");

        const RealT iq_endpoint  = std::sqrt(4.0 - 0.2 * 0.2);
        success                 *= initializationRejectedAtomically(makeData(), iq_endpoint, 0.2, 1.0, "reactive-current endpoint");
        success                 *= initializationRejectedAtomically(makeData(), -iq_endpoint, 0.2, 1.0, "negative reactive-current endpoint");

        auto q_priority                        = makeData();
        q_priority.parameters[Params::Pqflag]  = false;
        success                               *= initializationRejectedAtomically(q_priority, 0.2, iq_endpoint, 1.0, "Q-priority active-current endpoint");

        auto pord_limit                      = makeData();
        pord_limit.parameters[Params::Pmax]  = 0.25;
        success                             *= initializationRejectedAtomically(pord_limit, 0.0, 0.5, 1.0, "recovered PORD above Pmax");
        success                             *= initializationRejectedAtomically(makeData(), 0.0, 1.0e-6, 1.0, "recovered PORD below Pmin");

        auto q_endpoint                       = makeData();
        q_endpoint.parameters[Params::QFlag]  = true;
        q_endpoint.parameters[Params::VFlag]  = true;
        q_endpoint.parameters[Params::Kqi]    = 0.4;
        success                              *= initializationRejectedAtomically(q_endpoint, 0.0, 0.2, 1.0, "QGEN at Qmax", 0.2, 1.0);
        success                              *= initializationRejectedAtomically(q_endpoint, 0.0, 0.2, 1.0, "QGEN at Qmin", 0.2, -1.0);

        auto v_endpoint                     = q_endpoint;
        v_endpoint.parameters[Params::Kqi]  = 0.0;
        v_endpoint.parameters[Params::Kvi]  = 0.5;
        success                            *= initializationRejectedAtomically(v_endpoint, 0.0, 0.2, 1.2, "voltage at Vmax", 0.24, 0.0);

        // At Vmin the saturated Q-PI output equals the measured voltage, so
        // this boundary point is a consistent equilibrium and initializes.
        Fixture<ScalarT> v_boundary(v_endpoint, 0.8);
        v_boundary.attachAllInputs();
        v_boundary.input(Ext::PE)  = 0.16;
        success                   *= v_boundary.initialize(0.0, 0.2);
        success                   *= (v_boundary.evaluate() == 0);
        success                   *= allResidualsZero(v_boundary.reecb);

        auto zero_power                        = makeData();
        zero_power.parameters[Params::PfFlag]  = true;
        success                               *= initializationRejectedAtomically(zero_power, 0.1, 0.2, 1.0, "power-factor target at zero active power", 0.0, 0.1);

        auto pf_resolution                        = makeData();
        pf_resolution.parameters[Params::PfFlag]  = true;
        success                                  *= initializationRejectedAtomically(
            pf_resolution, 0.1, 0.2, 0.01, "unrepresentable power-factor reference", 1.0e-8);

        success *= initializationRejectedAtomically(makeData(), 0.0, 0.2, 0.0, "zero terminal voltage");
        success *= initializationRejectedAtomically(makeData(), 0.0, 0.2, 1.0, "nonfinite PE", std::numeric_limits<RealT>::infinity(), 0.0);

        auto collapsed                      = makeData();
        collapsed.parameters[Params::QFlag] = true;
        collapsed.parameters[Params::VFlag] = true;
        collapsed.parameters[Params::Kqi]   = 0.4;
        collapsed.parameters[Params::Kvi]   = 0.5;
        collapsed.parameters[Params::Qmin]  = 0.0;
        collapsed.parameters[Params::Qmax]  = 0.0;
        collapsed.parameters[Params::Vmin]  = 1.0;
        collapsed.parameters[Params::Vmax]  = 1.0;
        Fixture<ScalarT> collapsed_fixture(collapsed);
        collapsed_fixture.attachAllInputs();
        collapsed_fixture.input(Ext::PE)  = 0.2;
        success                          *= collapsed_fixture.initialize(0.0, 0.2);
        success                          *= (collapsed_fixture.evaluate() == 0);
        success                          *= allResidualsZero(collapsed_fixture.reecb);

        auto collapsed_q                      = collapsed;
        collapsed_q.parameters[Params::Vmin]  = 0.8;
        collapsed_q.parameters[Params::Vmax]  = 1.2;
        collapsed_q.parameters[Params::Kvi]   = 0.0;
        success                              *= initializationRejectedAtomically(collapsed_q, 0.0, 0.2, 1.0, "collapsed Q limit away from equilibrium", 0.2, 0.1);

        auto collapsed_v                      = collapsed;
        collapsed_v.parameters[Params::Qmin]  = -1.0;
        collapsed_v.parameters[Params::Qmax]  = 1.0;
        collapsed_v.parameters[Params::Kqi]   = 0.0;
        collapsed_v.parameters[Params::Vmin]  = 1.1;
        collapsed_v.parameters[Params::Vmax]  = 1.1;
        success                              *= initializationRejectedAtomically(collapsed_v, 0.0, 0.2, 1.0, "collapsed V limit away from equilibrium", 0.2, 0.0);

        auto zero_gains                      = makeData();
        zero_gains.parameters[Params::QFlag] = true;
        zero_gains.parameters[Params::VFlag] = true;
        zero_gains.parameters[Params::Kqi]   = 0.0;
        zero_gains.parameters[Params::Kvi]   = 0.0;
        Fixture<ScalarT> unconstrained(zero_gains);
        unconstrained.attachAllInputs();
        unconstrained.input(Ext::PE)    = 0.2;
        unconstrained.input(Ext::QGEN)  = 4.0;
        success                        *= unconstrained.initialize(0.0, 0.2);
        success                        *= (unconstrained.evaluate() == 0);
        success                        *= allResidualsZero(unconstrained.reecb);

        auto unattached_data                      = makeData();
        unattached_data.parameters[Params::QFlag] = true;
        unattached_data.parameters[Params::VFlag] = true;
        unattached_data.parameters[Params::Kqi]   = 0.4;
        unattached_data.parameters[Params::Kvi]   = 0.5;
        unattached_data.parameters[Params::kqv]   = 1.0;
        unattached_data.parameters.erase(Params::Vref0);
        Fixture<ScalarT> unattached(unattached_data, 1.2);
        success *= unattached.prepare(0.05, 0.2);
        setControlState(unattached.reecb);
        success                    *= (unattached.evaluate() == 0);
        const auto residual_before  = snapshot(unattached.reecb.getResidual());
        success                    *= (unattached.reecb.initialize() != 0);
        success                    *= (unattached.evaluate() == 0);
        success                    *= vectorUnchanged(unattached.reecb.getResidual(), residual_before, "unattached residual");

        return success.report(__func__);
      }

      /// Fixed near-endpoint commands check that the private smooth inverse reproduces each requested command without an artificial offset.
      TestOutcome initializationExactness()
      {
        TestStatus success = true;

        auto data                     = makeData();
        data.parameters[Params::Pmin] = -1.0;
        data.parameters[Params::Pmax] = 3.0;

        struct ExactnessCase
        {
          RealT       ipcmd;
          RealT       pord;
          const char* label;
        };

        const std::array<ExactnessCase, 3> cases{{
            {1.0e-6, -0.034728131800926182, "near lower active-current limit"},
            {0.2, 0.2, "interior active-current command"},
            {1.999999, 2.0347281318012689, "near upper active-current limit"},
        }};

        for (const auto& test_case : cases)
        {
          Fixture<ScalarT> fixture(data);
          success *= fixture.initialize(0.0, test_case.ipcmd);
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.y().getData()[index(Vars::PORD)], test_case.pord, test_case.label);
          success *= (fixture.ipcmd() == test_case.ipcmd);
          success *= allResidualsZero(fixture.reecb);
        }

        const std::array<ExactnessCase, 2> pord_boundaries{{
            {clampEdgeOffset(), 0.0, "PORD at Pmin"},
            {1.0, 1.0, "PORD at Pmax"},
        }};

        for (const auto& test_case : pord_boundaries)
        {
          Fixture<ScalarT> fixture(makeData());
          success *= fixture.initialize(0.0, test_case.ipcmd);
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.y().getData()[index(Vars::PORD)], test_case.pord, test_case.label);
          success *= allResidualsZero(fixture.reecb);
        }

        const RealT pmax                   = 0.2 - 0.5 * kTol;
        auto        near_pmax              = makeData();
        near_pmax.parameters[Params::Pmax] = pmax;
        Fixture<ScalarT> exact(near_pmax);
        success                    *= exact.initialize(0.0, 0.2);
        success                    *= (exact.evaluate() == 0);
        const RealT recovered_pord  = static_cast<RealT>(exact.reecb.y().getData()[index(Vars::PORD)]);
        success                    *= (recovered_pord > pmax);
        success                    *= (recovered_pord < pmax + kTol);
        success                    *= allResidualsZero(exact.reecb);

        const RealT      iqmax = std::sqrt(3.96);
        Fixture<ScalarT> reactive(data);
        success *= reactive.initialize(iqmax - 1.0e-6, 0.2);
        success *= (reactive.evaluate() == 0);
        success *= scalarMatches(reactive.reecb.y().getData()[index(Vars::QV)], 2.0247030060145095, "near upper reactive-current limit");
        success *= allResidualsZero(reactive.reecb);

        return success.report(__func__);
      }

      /// One independently calculated literal answer key checks all ten residual rows at a rich non-equilibrium state.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeDynamicData(), 0.9, 0.4);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.25, 0.4);
        setAnswerKeyState(fixture.reecb);
        success *= (fixture.evaluate() == 0);

        const std::array<ExpectedResidual, 10> expected{{
            {Vars::VMEAS, "VMEAS", 0.24000000000000021},  // -VMEAS' + (VT - VMEAS) / Trv
            {Vars::PMEAS, "PMEAS", 0.14499999999999982},  // -PMEAS' + (kbase PE - PMEAS) / Tp
            {Vars::XPIQ, "XPIQ", 0.083249747972647115},   // -XPIQ' + sQPI sdip antiwindup(Kqp eq + XPIQ, Kqi eq; Vmin, Vmax)
            {Vars::XPIV, "XPIV", -0.095000000000000001},  // -XPIV' + sQ sdip antiwindup(Kvp epiv + XPIV, Kvi epiv; -Iqmax, Iqmax)
            {Vars::QV, "QV", -0.050000000000000003},      // -QV' + sQoff sdip (qref / vsafe - QV) / Tiq
            {Vars::PORD, "PORD", 0.25999999999999973},    // -PORD' + sdip antiwindup(PORD, rpord; Pmin, Pmax)
            {Vars::VT, "VT", -0.029999999999999916},      // -VT^2 + Vr^2 + Vi^2
            {Vars::ILMAX, "ILMAX", 0.16999999999999993},  // -ILMAX |ILMAX| + Imax^2 - sPQ (kbase IPCMD)^2 - sPQoff (kbase IQCMD)^2
            {Vars::IQCMD, "IQCMD", -0.76999943561644268}, // -kbase IQCMD + clamp(iqraw; -Iqmax, Iqmax)
            {Vars::IPCMD, "IPCMD", -0.11578947368421055}, // -kbase IPCMD + clamp(PORD / vsafe; 0, Ipmax)
        }};

        const auto* residual = fixture.reecb.getResidual().getData();
        for (const auto& answer : expected)
        {
          success *= scalarMatches(residual[index(answer.row)], answer.value, answer.name);
        }

        return success.report(__func__);
      }

      /// Every valid selector combination initializes attached and unattached signals to a
      /// zero-residual state; power-factor control with the direct-voltage reference is rejected.
      TestOutcome selectorConfigurations()
      {
        TestStatus success = true;

        noteExpectedLogs("Testing REECB selector configurations. "
                         "Rejection of the power-factor direct-voltage combinations is expected.");

        for (const bool pf : {false, true})
        {
          for (const bool voltage : {false, true})
          {
            for (const bool reactive : {false, true})
            {
              for (const bool p_priority : {false, true})
              {
                auto data                       = makeData();
                data.parameters[Params::PfFlag] = pf;
                data.parameters[Params::VFlag]  = voltage;
                data.parameters[Params::QFlag]  = reactive;
                data.parameters[Params::Pqflag] = p_priority;
                data.parameters[Params::Kqi]    = reactive && voltage ? 0.4 : 0.0;
                data.parameters[Params::Kvi]    = reactive ? 0.5 : 0.0;

                for (const bool attached : {false, true})
                {
                  Fixture<ScalarT> fixture(data);
                  if (attached)
                  {
                    fixture.attachAllInputs(7.0);
                    fixture.input(Ext::PE)   = 0.2;
                    fixture.input(Ext::QGEN) = 0.1;
                  }

                  if (pf && reactive && !voltage)
                  {
                    success *= (fixture.reecb.verify() > 0);
                    success *= !fixture.initialize(0.1, 0.2);
                    continue;
                  }

                  success *= fixture.initialize(0.1, 0.2);
                  success *= (fixture.evaluate() == 0);
                  success *= allResidualsZero(fixture.reecb);
                  success *= (fixture.iqcmd() == 0.1);
                  success *= (fixture.ipcmd() == 0.2);

                  const auto* y               = fixture.reecb.y().getData();
                  const RealT expected_ilmax  = p_priority ? std::sqrt(3.96) : std::sqrt(3.99);
                  success                    *= scalarMatches(y[index(Vars::ILMAX)], expected_ilmax, "selector ILMAX");

                  if (reactive)
                  {
                    success *= std::abs(y[index(Vars::XPIV)]) > kTol;
                    if (voltage)
                    {
                      success *= std::abs(y[index(Vars::XPIQ)]) > kTol;
                    }
                  }
                  else
                  {
                    success *= std::abs(y[index(Vars::QV)]) > kTol;
                  }

                  if (attached)
                  {
                    success                     *= (fixture.input(Ext::PE) == 0.2);
                    success                     *= (fixture.input(Ext::QGEN) == 0.1);
                    const RealT qref             = reactive && !voltage ? 1.0 : 0.1;
                    const RealT expected_pfaref  = pf ? 0.4636476090008061 : 0.0;
                    success                     *= scalarMatches(fixture.input(Ext::QEXT), qref, "selector QEXT publication");
                    success                     *= scalarMatches(fixture.input(Ext::PFAREF), expected_pfaref, "selector PFAREF publication");
                    success                     *= scalarMatches(fixture.input(Ext::PREF), 0.2, "selector PREF publication");
                  }
                }
              }
            }
          }
        }

        return success.report(__func__);
      }

      /// Direct-voltage mode consumes and publishes the Volt/VAr reference without
      /// power-base conversion; reactive modes convert on the nonidentity base.
      TestOutcome voltVarReferenceBase()
      {
        TestStatus success = true;

        {
          // 50 MVA component on the 100 MVA system selects direct-voltage mode.
          auto data                      = makeData();
          data.parameters[Params::mva]   = 50.0;
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = false;
          data.parameters[Params::Kvi]   = 0.5;

          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          success *= fixture.initialize(0.1, 0.2);
          success *= scalarMatches(fixture.input(Ext::QEXT), 1.0, "published voltage reference");
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.reecb);

          // A raised external voltage reference enters the V-PI rate raw.
          fixture.input(Ext::QEXT)  = 1.02;
          success                  *= (fixture.evaluate() == 0);
          success                  *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIV)], 0.01, "unconverted voltage-reference rate");
        }

        {
          // The reactive-current lag keeps the power-base conversion.
          auto data                    = makeData();
          data.parameters[Params::mva] = 50.0;

          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          success *= fixture.initialize(0.1, 0.2);
          success *= scalarMatches(fixture.input(Ext::QEXT), 0.1, "published system-base reactive power");
          success *= (fixture.evaluate() == 0);
          success *= allResidualsZero(fixture.reecb);

          fixture.input(Ext::QEXT)  = 0.11;
          success                  *= (fixture.evaluate() == 0);
          success                  *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::QV)], 1.0, "converted reactive-reference rate");
        }

        return success.report(__func__);
      }

      /// Literal selector, voltage-gate, limiter, and anti-windup cases check the reactive-control paths.
      TestOutcome reactiveControl()
      {
        TestStatus success = true;

        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = false;
          data.parameters[Params::kqv]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 0.4;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::QV)] = 0.1;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::QV)], 2.3333333333333335, "constant-Q lag");
        }

        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = true;
          data.parameters[Params::kqv]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 0.1;
          fixture.input(Ext::QGEN)  = -0.05;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::XPIQ)] = 1.0;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIQ)], 0.11999999999996271, "Q-control integral rate");
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::QV)], 0.0, "bypassed Q lag");
        }

        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = false;
          data.parameters[Params::kqv]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 1.05;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIQ)], 0.0, "bypassed Q-control integrator");
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIV)], 0.025, "voltage-control integral rate");
        }

        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = false;
          data.parameters[Params::kqv]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 0.4;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          auto* y            = fixture.reecb.y().getData();
          y[index(Vars::QV)] = 0.1;
          y[index(Vars::VT)] = 0.5;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::QV)], 0.0, "constant-Q voltage-band gate");
        }

        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = true;
          data.parameters[Params::kqv]    = 0.0;
          data.parameters[Params::Qmin]   = -2.0;
          data.parameters[Params::Qmax]   = 2.0;
          data.parameters[Params::Kqp]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 0.5;
          fixture.input(Ext::QGEN)  = 0.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          auto* y              = fixture.reecb.y().getData();
          y[index(Vars::XPIQ)] = 1.0;
          y[index(Vars::VT)]   = 0.5;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIQ)], 0.0, "Q-integrator voltage-band gate");
        }

        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = false;
          data.parameters[Params::kqv]    = 0.0;
          data.parameters[Params::Kvp]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 1.1;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          auto* y              = fixture.reecb.y().getData();
          y[index(Vars::XPIV)] = 0.0;
          y[index(Vars::VT)]   = 0.5;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIV)], 0.0, "V-integrator voltage-band gate");
        }

        struct ReactiveCase
        {
          RealT       input;
          RealT       state;
          RealT       expected;
          const char* label;
        };

        const std::array<ReactiveCase, 5> q_limit_cases{{
            {-1.0, 1.0, -0.27999999999999997, "Q reference below Qmin"},
            {-0.7, 1.0, 0.4 * (-0.7 + clampEdgeOffset()), "Q reference at Qmin"},
            {0.1, 1.0, 0.039999999999999994, "Q reference inside limits"},
            {0.8, 1.0, 0.4 * (0.8 - clampEdgeOffset()), "Q reference at Qmax"},
            {1.0, 1.0, 0.32000000000000006, "Q reference above Qmax"},
        }};

        for (const auto& test_case : q_limit_cases)
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = true;
          data.parameters[Params::kqv]    = 0.0;
          data.parameters[Params::Kqp]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = test_case.input / 2.0;
          fixture.input(Ext::QGEN)  = 0.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::XPIQ)] = test_case.state;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIQ)], test_case.expected, test_case.label);
        }

        const std::array<ReactiveCase, 6> q_windup_cases{{
            {1.0, 2.0, 0.0, "outward Q-integrator rate above Vmax"},
            {1.0, 1.3, 0.2, "outward Q-integrator rate at Vmax"},
            {-1.0, 2.0, -0.4, "restoring Q-integrator rate above Vmax"},
            {-1.0, 0.0, 0.0, "outward Q-integrator rate below Vmin"},
            {-1.0, 0.7, -0.2, "outward Q-integrator rate at Vmin"},
            {1.0, 0.0, 0.4, "restoring Q-integrator rate below Vmin"},
        }};

        for (const auto& test_case : q_windup_cases)
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = true;
          data.parameters[Params::kqv]    = 0.0;
          data.parameters[Params::Qmin]   = -2.0;
          data.parameters[Params::Qmax]   = 2.0;
          data.parameters[Params::Kqp]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = test_case.input / 2.0;
          fixture.input(Ext::QGEN)  = 0.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::XPIQ)] = test_case.state;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIQ)], test_case.expected, test_case.label);
        }

        const std::array<ReactiveCase, 5> v_limit_cases{{
            {0.0, 0.0, -0.30000000000000004, "V-PI input below Vmin"},
            {0.7, 0.0, -0.3 + clampEdgeOffset(), "V-PI input at Vmin"},
            {1.0, 0.0, 0.0, "V-PI input inside limits"},
            {1.3, 0.0, 0.3 - clampEdgeOffset(), "V-PI input at Vmax"},
            {2.0, 0.0, 0.30000000000000004, "V-PI input above Vmax"},
        }};

        for (const auto& test_case : v_limit_cases)
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = true;
          data.parameters[Params::kqv]    = 0.0;
          data.parameters[Params::Kqp]    = 0.0;
          data.parameters[Params::Kvi]    = 1.0;
          data.parameters[Params::Kvp]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = 0.0;
          fixture.input(Ext::QGEN)  = 0.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          auto* y              = fixture.reecb.y().getData();
          y[index(Vars::XPIQ)] = test_case.input;
          y[index(Vars::XPIV)] = test_case.state;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIV)], test_case.expected, test_case.label);
        }

        const std::array<ReactiveCase, 6> v_windup_cases{{
            {0.5, 2.0, 0.0, "outward V-integrator rate above Iqmax"},
            {0.5, 1.4, 0.25, "outward V-integrator rate at Iqmax"},
            {-0.5, 2.0, -0.5, "restoring V-integrator rate above Iqmax"},
            {-0.5, -2.0, 0.0, "outward V-integrator rate below negative Iqmax"},
            {-0.5, -1.4, -0.25, "outward V-integrator rate at negative Iqmax"},
            {0.5, -2.0, 0.5, "restoring V-integrator rate below negative Iqmax"},
        }};

        for (const auto& test_case : v_windup_cases)
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::PfFlag] = false;
          data.parameters[Params::QFlag]  = true;
          data.parameters[Params::VFlag]  = false;
          data.parameters[Params::kqv]    = 0.0;
          data.parameters[Params::Kvp]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::QEXT)  = test_case.input > 0.0 ? 2.0 : 0.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::XPIV)] = test_case.state;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::XPIV)], test_case.expected, test_case.label);
        }

        // With QFlag = 0 and a zero lag state the reactive command is the
        // injection alone, so the IQCMD row reads iqv = clamp(kqv deadband2(
        // Vref0 - VMEAS; dbd1, dbd2); Iql1, Iqh1) directly. The wide deadband
        // and limit gaps keep every smoothing tail below kTol except at the
        // tested edges, which sit exactly one clampEdgeOffset() inside.
        const std::array<ReactiveCase, 7> injection_cases{{
            {1.75, 0.0, -0.4, "injection saturated at Iql1"},
            {1.2, 0.0, -clampEdgeOffset(), "injection at lower deadband breakpoint"},
            {1.0, 0.0, 0.0, "injection inside deadband"},
            {0.75, 0.0, clampEdgeOffset(), "injection at upper deadband breakpoint"},
            {0.55, 0.0, 0.2, "injection passthrough above deadband"},
            {0.25, 0.0, 0.5 - clampEdgeOffset(), "injection at Iqh1 edge"},
            {0.1, 0.0, 0.5, "injection saturated at Iqh1"},
        }};

        for (const auto& test_case : injection_cases)
        {
          auto data                     = makeData();
          data.parameters[Params::kqv]  = 1.0;
          data.parameters[Params::dbd1] = -0.2;
          data.parameters[Params::dbd2] = 0.25;
          data.parameters[Params::Iql1] = -0.4;
          data.parameters[Params::Iqh1] = 0.5;
          Fixture<ScalarT> fixture(data);
          success *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::VMEAS)] = test_case.input;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::IQCMD)], test_case.expected, test_case.label);
        }

        return success.report(__func__);
      }

      /// Literal ramp, voltage-gate, anti-windup, command-limit, current-circle, and signed-continuation cases check the active/current paths.
      TestOutcome activeCurrentControl()
      {
        TestStatus success = true;

        struct ScalarCase
        {
          RealT       input;
          RealT       expected;
          const char* label;
        };

        const std::array<ScalarCase, 5> rate_cases{{
            {-1.0, -0.5, "lower PORD ramp saturation"},
            {-0.5, -0.5 + clampEdgeOffset(), "lower PORD ramp boundary"},
            {0.2, 0.19999999999999996, "interior PORD rate"},
            {0.6, 0.6 - clampEdgeOffset(), "upper PORD ramp boundary"},
            {1.0, 0.6, "upper PORD ramp saturation"},
        }};

        for (const auto& test_case : rate_cases)
        {
          Fixture<ScalarT> fixture(makeDynamicData());
          fixture.attachAllInputs();
          fixture.input(Ext::PREF)  = (0.65 + 0.25 * test_case.input) / 2.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::PORD)] = 0.65;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::PORD)], test_case.expected, test_case.label);
        }

        const std::array<ScalarCase, 5> gate_cases{{
            {0.5, 0.0, "PORD below voltage band"},
            {0.7, 0.1, "PORD at lower voltage threshold"},
            {1.0, 0.2, "PORD inside voltage band"},
            {1.2, 0.1, "PORD at upper voltage threshold"},
            {1.4, 0.0, "PORD above voltage band"},
        }};

        for (const auto& test_case : gate_cases)
        {
          Fixture<ScalarT> fixture(makeDynamicData());
          fixture.attachAllInputs();
          fixture.input(Ext::PREF)  = 0.35;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          auto* y              = fixture.reecb.y().getData();
          y[index(Vars::PORD)] = 0.65;
          y[index(Vars::VT)]   = test_case.input;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::PORD)], test_case.expected, test_case.label);
        }

        struct WindupCase
        {
          RealT       pord;
          RealT       raw_rate;
          RealT       expected;
          const char* label;
        };

        const std::array<WindupCase, 6> windup_cases{{
            {2.0, 1.0, 0.0, "outward rate above Pmax"},
            {1.4, 1.0, 0.3, "outward rate at Pmax"},
            {2.0, -1.0, -0.5, "restoring rate above Pmax"},
            {-1.0, -1.0, 0.0, "outward rate below Pmin"},
            {0.1, -1.0, -0.25, "outward rate at Pmin"},
            {-1.0, 1.0, 0.6, "restoring rate below Pmin"},
        }};

        for (const auto& test_case : windup_cases)
        {
          Fixture<ScalarT> fixture(makeDynamicData());
          fixture.attachAllInputs();
          fixture.input(Ext::PREF)  = (test_case.pord + 0.25 * test_case.raw_rate) / 2.0;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          fixture.reecb.y().getData()[index(Vars::PORD)] = test_case.pord;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::PORD)], test_case.expected, test_case.label);
        }

        const std::array<ScalarCase, 5> iq_limit_cases{{
            {-1.0, -0.4, "IQCMD below low-priority limit"},
            {-0.4, -0.4 + clampEdgeOffset(), "IQCMD at negative low-priority limit"},
            {0.0, 0.0, "IQCMD inside low-priority limits"},
            {0.4, 0.4 - clampEdgeOffset(), "IQCMD at positive low-priority limit"},
            {1.0, 0.4, "IQCMD above low-priority limit"},
        }};

        for (const auto& test_case : iq_limit_cases)
        {
          auto data                       = makeData();
          data.parameters[Params::QFlag]  = false;
          data.parameters[Params::Pqflag] = true;
          Fixture<ScalarT> fixture(data);
          success *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          auto* y               = fixture.reecb.y().getData();
          y[index(Vars::ILMAX)] = 0.4;
          y[index(Vars::IQCMD)] = 0.0;
          y[index(Vars::QV)]    = test_case.input;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::IQCMD)], test_case.expected, test_case.label);
        }

        const std::array<ScalarCase, 5> ip_limit_cases{{
            {-1.0, 0.0, "IPCMD below zero"},
            {0.0, clampEdgeOffset(), "IPCMD at zero"},
            {0.2, 0.2, "IPCMD inside low-priority limits"},
            {0.4, 0.4 - clampEdgeOffset(), "IPCMD at low-priority limit"},
            {1.0, 0.4, "IPCMD above low-priority limit"},
        }};

        for (const auto& test_case : ip_limit_cases)
        {
          auto data                       = makeData();
          data.parameters[Params::Pqflag] = false;
          Fixture<ScalarT> fixture(data);
          success *= fixture.prepare(0.2, 0.0);
          setControlState(fixture.reecb);
          auto* y               = fixture.reecb.y().getData();
          y[index(Vars::ILMAX)] = 0.4;
          y[index(Vars::IPCMD)] = 0.0;
          y[index(Vars::PORD)]  = test_case.input;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::IPCMD)], test_case.expected, test_case.label);
        }

        for (const bool p_priority : {false, true})
        {
          auto data                       = makeDynamicData();
          data.parameters[Params::Pqflag] = p_priority;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          success *= fixture.prepare(0.25, 0.4);
          setAnswerKeyState(fixture.reecb);
          success              *= (fixture.evaluate() == 0);
          const RealT expected  = p_priority ? 0.16999999999999993 : 0.56;
          success              *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::ILMAX)], expected, "priority-circle residual");
        }

        {
          auto data                     = makeData();
          data.parameters[Params::Imax] = 1.0;
          Fixture<ScalarT> fixture(data);
          const RealT      ipcmd  = std::sqrt(1.0 - 1.0e-12);
          success                *= fixture.prepare(0.0, ipcmd);
          setControlState(fixture.reecb);
          auto* y               = fixture.reecb.y().getData();
          y[index(Vars::ILMAX)] = 1.0e-6;
          y[index(Vars::IPCMD)] = ipcmd;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::ILMAX)], 0.0, "near-zero ILMAX residual");
        }

        for (const bool p_priority : {false, true})
        {
          auto data                       = makeData();
          data.parameters[Params::Pqflag] = p_priority;
          data.parameters[Params::Imax]   = 1.0;

          const RealT      iqcmd = p_priority ? 0.2 : 1.1;
          const RealT      ipcmd = p_priority ? 1.1 : 0.2;
          Fixture<ScalarT> fixture(data);
          success *= fixture.prepare(iqcmd, ipcmd);
          setControlState(fixture.reecb);

          auto* y                                             = fixture.reecb.y().getData();
          y[index(Vars::ILMAX)]                               = -std::sqrt(0.21);
          y[p_priority ? index(Vars::QV) : index(Vars::PORD)] = 0.3;
          fixture.reecb.y().setDataUpdated();
          success                          *= (fixture.evaluate() == 0);
          success                          *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::ILMAX)], 0.0, "negative ILMAX continuation");
          const auto  low_priority_row      = p_priority ? Vars::IQCMD : Vars::IPCMD;
          const RealT negative_limit_value  = fixture.reecb.getResidual().getData()[index(low_priority_row)];
          success                          *= scalarMatches(negative_limit_value, 0.1, "negative ILMAX command bound");
          for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
          {
            success *= std::isfinite(fixture.reecb.getResidual().getData()[row]);
          }

          y[index(Vars::ILMAX)] = std::sqrt(0.21);
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(low_priority_row)], negative_limit_value, "signed ILMAX bound parity");

          y[index(Vars::ILMAX)] = 0.0;
          fixture.reecb.y().setDataUpdated();
          success *= (fixture.evaluate() == 0);
          success *= scalarMatches(fixture.reecb.getResidual().getData()[index(Vars::ILMAX)], -0.21, "zero ILMAX continuation");
          for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
          {
            success *= std::isfinite(fixture.reecb.getResidual().getData()[row]);
          }
        }

        return success.report(__func__);
      }

      /// Fixed positive and negative DependencyTracking coefficients provide the oracle before each configuration is compared with Enzyme.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        constexpr RealT alpha = 0.7;
        for (const bool pf : {false, true})
        {
          for (const bool voltage : {false, true})
          {
            for (const bool reactive : {false, true})
            {
              for (const bool p_priority : {false, true})
              {
                if (pf && reactive && !voltage)
                {
                  continue;
                }

                auto data                       = makeJacobianData();
                data.parameters[Params::PfFlag] = pf;
                data.parameters[Params::VFlag]  = voltage;
                data.parameters[Params::QFlag]  = reactive;
                data.parameters[Params::Pqflag] = p_priority;

                const auto dependency    = dependencyTrackingJacobian(data, alpha, success);
                const auto bus_vr        = index(Vars::MAXIMUM);
                const auto bus_vi        = bus_vr + 1;
                const auto pe_column     = bus_vi + 1 + index(Ext::PE);
                const auto qgen_column   = bus_vi + 1 + index(Ext::QGEN);
                const auto qext_column   = bus_vi + 1 + index(Ext::QEXT);
                const auto pfaref_column = bus_vi + 1 + index(Ext::PFAREF);
                const auto pref_column   = bus_vi + 1 + index(Ext::PREF);

                success *= derivativeMatches(dependency, Vars::VMEAS, Vars::VMEAS, -5.7, "VMEAS diagonal");
                success *= derivativeMatches(dependency, Vars::VMEAS, Vars::VT, 5.0, "VMEAS-VT");
                success *= derivativeMatches(dependency, Vars::PMEAS, Vars::PMEAS, -3.2, "PMEAS diagonal");
                success *= derivativeMatches(dependency, Vars::PMEAS, pe_column, 5.0, "PMEAS-PE");
                success *= derivativeMatches(dependency, Vars::XPIQ, Vars::XPIQ, -alpha, "XPIQ diagonal");
                success *= derivativeMatches(dependency, Vars::XPIV, Vars::XPIV, -alpha, "XPIV diagonal");
                success *= derivativeMatches(dependency, Vars::QV, Vars::QV, -alpha - (reactive ? 0.0 : 1.0 / 0.3), "QV diagonal"); // Tiq = 0.3
                success *= derivativeMatches(dependency, Vars::PORD, Vars::PORD, -4.7, "PORD diagonal");
                success *= derivativeMatches(dependency, Vars::PORD, pref_column, 8.0, "PORD-PREF");
                success *= derivativeMatches(dependency, Vars::VT, Vars::VT, -2.0, "VT diagonal");
                success *= derivativeMatches(dependency, Vars::VT, bus_vr, 1.8, "VT-Vr");
                success *= derivativeMatches(dependency, Vars::VT, bus_vi, 0.8, "VT-Vi");
                success *= derivativeMatches(dependency, Vars::ILMAX, Vars::ILMAX, -2.4, "ILMAX diagonal");
                success *= derivativeMatches(dependency, Vars::IQCMD, Vars::IQCMD, -2.0, "IQCMD diagonal");
                success *= derivativeMatches(dependency, Vars::IPCMD, Vars::IPCMD, -2.0, "IPCMD diagonal");
                success *= derivativeMatches(dependency, Vars::IPCMD, Vars::PORD, 1.0, "IPCMD-PORD");
                success *= derivativeMatches(dependency, Vars::IPCMD, Vars::VMEAS, -0.5, "IPCMD-VMEAS");
                success *= derivativeMatches(dependency, Vars::IQCMD, Vars::XPIV, reactive ? 1.0 : 0.0, "IQCMD-XPIV selector path");
                success *= derivativeMatches(dependency, Vars::IQCMD, Vars::QV, reactive ? 0.0 : 1.0, "IQCMD-QV selector path");
                success *= derivativeMatches(dependency, Vars::XPIQ, qgen_column, reactive && voltage ? -0.8 : 0.0, "XPIQ-QGEN selector path");

                // The direct-voltage coefficient carries no power-base factor,
                // while the cascaded path converts the reference to component base.
                RealT xpiv_qext = 0.0;
                if (reactive)
                {
                  xpiv_qext = voltage ? (pf ? 0.0 : 0.6) : 0.5;
                }
                success *= derivativeMatches(dependency, Vars::XPIV, qext_column, xpiv_qext, "XPIV-QEXT selector path");
                success *= derivativeMatches(dependency, Vars::QV, qext_column, !reactive && !pf ? 20.0 / 3.0 : 0.0, "QV-QEXT selector path");
                success *= derivativeMatches(dependency, Vars::QV, pfaref_column, !reactive && pf ? 25.0 / 3.0 : 0.0, "QV-PFAREF selector path");

                if (p_priority)
                {
                  success *= derivativeMatches(dependency, Vars::ILMAX, Vars::IPCMD, -1.6, "P-priority current-circle column");
                  success *= derivativeMatches(dependency, Vars::ILMAX, Vars::IQCMD, 0.0, "P-priority absent current-circle column");
                }
                else
                {
                  success *= derivativeMatches(dependency, Vars::ILMAX, Vars::IQCMD, -0.8, "Q-priority current-circle column");
                  success *= derivativeMatches(dependency, Vars::ILMAX, Vars::IPCMD, 0.0, "Q-priority absent current-circle column");
                }

#ifdef GRIDKIT_ENABLE_ENZYME
                const auto enzyme  = enzymeJacobian(data, alpha, success);
                success           *= jacobiansMatch(dependency, enzyme, index(Vars::MAXIMUM) + 2 + index(Ext::MAXIMUM));
#endif
              }
            }
          }
        }

        for (const bool p_priority : {false, true})
        {
          auto data                       = makeJacobianData();
          data.parameters[Params::Pqflag] = p_priority;

          const auto dependency  = dependencyTrackingJacobian(data, alpha, success, -1.2);
          success               *= derivativeMatches(dependency, Vars::ILMAX, Vars::ILMAX, -2.4, "negative ILMAX continuation");
#ifdef GRIDKIT_ENABLE_ENZYME
          const auto enzyme  = enzymeJacobian(data, alpha, success, -1.2);
          success           *= jacobiansMatch(dependency, enzyme, index(Vars::MAXIMUM) + 2 + index(Ext::MAXIMUM));
#endif
        }

        // The selector sweep zeroes kqv because the answer-key deadband tails
        // sit above kTol there. This configuration exercises the injection
        // derivative on its own: with QFlag = 0 the only VMEAS dependence of
        // the IQCMD row is iqv, evaluated on the deadband passthrough side
        // strictly inside the injection limits, where the chain collapses to
        // d(iqv)/d(VMEAS) = -kqv.
        {
          auto data                      = makeJacobianData();
          data.parameters[Params::QFlag] = false;
          data.parameters[Params::kqv]   = 1.0;
          data.parameters[Params::dbd1]  = -0.2;
          data.parameters[Params::dbd2]  = 0.25;
          data.parameters[Params::Vref0] = 1.5;

          const auto dependency  = dependencyTrackingJacobian(data, alpha, success);
          success               *= derivativeMatches(dependency, Vars::IQCMD, Vars::VMEAS, -1.0, "IQCMD-VMEAS injection path");
          success               *= derivativeMatches(dependency, Vars::IQCMD, Vars::QV, 1.0, "IQCMD-QV alongside injection");
#ifdef GRIDKIT_ENABLE_ENZYME
          const auto enzyme  = enzymeJacobian(data, alpha, success);
          success           *= jacobiansMatch(dependency, enzyme, index(Vars::MAXIMUM) + 2 + index(Ext::MAXIMUM));
#endif
        }

        return success.report(__func__);
      }

    private:
      using Params        = PhasorDynamics::Converter::ReecbParameters;
      using Vars          = PhasorDynamics::Converter::ReecbInternalVariables;
      using Ext           = PhasorDynamics::Converter::ReecbExternalVariables;
      using Mon           = PhasorDynamics::Converter::ReecbMonitorableVariables;
      using Data          = PhasorDynamics::Converter::ReecbData<RealT, IdxT>;
      using ReecbT        = PhasorDynamics::Converter::Reecb<ScalarT, IdxT>;
      using DependencyMap = DependencyTracking::Variable::DependencyMap;

      struct ExpectedResidual
      {
        Vars        row;
        const char* name;
        RealT       value;
      };

      static constexpr size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(Ext variable)
      {
        return static_cast<size_t>(variable);
      }

      template <typename T>
      class Fixture
      {
      private:
        std::array<T, index(Ext::MAXIMUM)>                                   input_values_{};
        std::array<IdxT, index(Ext::MAXIMUM)>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, index(Ext::MAXIMUM)> input_nodes_{};
        PhasorDynamics::SignalNode<T, IdxT>                                  iqcmd_node_;
        PhasorDynamics::SignalNode<T, IdxT>                                  ipcmd_node_;
        bool                                                                 commands_assigned_{true};

      public:
        explicit Fixture(const Data& data, RealT vr = 1.0, RealT vi = 0.0, RealT system_va_base = 100.0e6, bool assign_commands = true)
          : commands_assigned_(assign_commands), bus(static_cast<T>(vr), static_cast<T>(vi)), reecb(&bus, data)
        {
          reecb.setSystemBase(60.0, system_va_base);
          if (commands_assigned_)
          {
            reecb.getSignals().template assignSignalNode<Vars::IQCMD>(&iqcmd_node_);
            reecb.getSignals().template assignSignalNode<Vars::IPCMD>(&ipcmd_node_);
          }
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachAllInputs(RealT value = 0.0)
        {
          for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
          {
            input_values_[port]  = static_cast<T>(value);
            input_indices_[port] = reecb.size() + bus.size() + static_cast<IdxT>(port);
            input_nodes_[port].set(&input_values_[port], &input_indices_[port]);
          }

          auto& signals = reecb.getSignals();
          signals.template attachSignalNode<Ext::PE>(&input_nodes_[index(Ext::PE)]);
          signals.template attachSignalNode<Ext::QGEN>(&input_nodes_[index(Ext::QGEN)]);
          signals.template attachSignalNode<Ext::QEXT>(&input_nodes_[index(Ext::QEXT)]);
          signals.template attachSignalNode<Ext::PFAREF>(&input_nodes_[index(Ext::PFAREF)]);
          signals.template attachSignalNode<Ext::PREF>(&input_nodes_[index(Ext::PREF)]);
        }

        void setCommands(RealT iqcmd, RealT ipcmd)
        {
          if (commands_assigned_)
          {
            iqcmd_node_.init(static_cast<T>(iqcmd));
            ipcmd_node_.init(static_cast<T>(ipcmd));
          }
          else
          {
            auto* y               = reecb.y().getData();
            y[index(Vars::IQCMD)] = static_cast<T>(iqcmd);
            y[index(Vars::IPCMD)] = static_cast<T>(ipcmd);
            reecb.y().setDataUpdated();
          }
        }

        bool prepare(RealT iqcmd, RealT ipcmd)
        {
          const bool ready = (bus.allocate() == 0) && (reecb.allocate() == 0)
                             && (reecb.verify() == 0) && (bus.initialize() == 0);
          if (!ready)
          {
            std::cout << "REECB fixture preparation failed\n";
            return false;
          }
          setCommands(iqcmd, ipcmd);
          return true;
        }

        bool initialize(RealT iqcmd, RealT ipcmd)
        {
          if (!prepare(iqcmd, ipcmd))
          {
            return false;
          }
          if (reecb.initialize() == 0)
          {
            return true;
          }
          std::cout << "REECB initialization failed\n";
          return false;
        }

        int evaluate()
        {
          return reecb.evaluateResidual();
        }

        T iqcmd() const
        {
          return commands_assigned_ ? iqcmd_node_.read() : reecb.y().getData()[index(Vars::IQCMD)];
        }

        T ipcmd() const
        {
          return commands_assigned_ ? ipcmd_node_.read() : reecb.y().getData()[index(Vars::IPCMD)];
        }

        T& input(Ext variable)
        {
          return input_values_[index(variable)];
        }

        IdxT inputIndex(Ext variable) const
        {
          return input_indices_[index(variable)];
        }

        PhasorDynamics::Bus<T, IdxT>              bus;
        PhasorDynamics::Converter::Reecb<T, IdxT> reecb;
      };

      Data makeData() const
      {
        Data data;
        data.device_class          = "Reecb";
        data.disambiguation_string = "reecb_test";
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);

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

      Data makeJacobianData() const
      {
        auto data                      = makeDynamicData();
        data.parameters[Params::Vref0] = std::sqrt(0.97);
        data.parameters[Params::kqv]   = 0.0;
        data.parameters[Params::Qmin]  = -2.0;
        data.parameters[Params::Qmax]  = 2.0;
        return data;
      }

      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::PE)     = 0.3;
        fixture.input(Ext::QGEN)   = -0.1;
        fixture.input(Ext::QEXT)   = 0.2;
        fixture.input(Ext::PFAREF) = 0.15;
        fixture.input(Ext::PREF)   = 0.35;
      }

      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Converter::Reecb<T, IdxT>& reecb) const
      {
        auto* y               = reecb.y().getData();
        y[index(Vars::VMEAS)] = 0.95;
        y[index(Vars::PMEAS)] = 0.55;
        y[index(Vars::XPIQ)]  = 0.10;
        y[index(Vars::XPIV)]  = -0.05;
        y[index(Vars::QV)]    = 0.30;
        y[index(Vars::PORD)]  = 0.65;
        y[index(Vars::VT)]    = 1.00;
        y[index(Vars::ILMAX)] = 1.20;
        y[index(Vars::IQCMD)] = 0.25;
        y[index(Vars::IPCMD)] = 0.40;

        auto* yp               = reecb.yp().getData();
        yp[index(Vars::VMEAS)] = 0.01;
        yp[index(Vars::PMEAS)] = -0.02;
        yp[index(Vars::XPIQ)]  = 0.03;
        yp[index(Vars::XPIV)]  = -0.03;
        yp[index(Vars::QV)]    = 0.05;
        yp[index(Vars::PORD)]  = -0.06;
        reecb.y().setDataUpdated();
        reecb.yp().setDataUpdated();
      }

      template <typename T>
      void setControlState(PhasorDynamics::Converter::Reecb<T, IdxT>& reecb) const
      {
        auto* y               = reecb.y().getData();
        y[index(Vars::VMEAS)] = 1.0;
        y[index(Vars::PMEAS)] = 0.0;
        y[index(Vars::XPIQ)]  = 0.0;
        y[index(Vars::XPIV)]  = 0.0;
        y[index(Vars::QV)]    = 0.0;
        y[index(Vars::PORD)]  = 0.5;
        y[index(Vars::VT)]    = 1.0;
        y[index(Vars::ILMAX)] = 1.4;
        reecb.yp().setToConst(static_cast<T>(0.0));
        reecb.y().setDataUpdated();
      }

      template <typename ValueT>
      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus, Params parameter, ValueT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> model(&bus, data);
        model.setSystemBase(60.0, 100.0e6);
        return model.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected(PhasorDynamics::Bus<ScalarT, IdxT>& bus) const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT>       node;
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> model(&bus, makeData());
        model.setSystemBase(60.0, 100.0e6);
        model.getSignals().template attachSignalNode<variable>(&node);
        return model.verify() > 0;
      }

      bool initializationRejectedAtomically(const Data& data, RealT iqcmd, RealT ipcmd, RealT voltage, const char* label, RealT pe = std::numeric_limits<RealT>::quiet_NaN(), RealT qgen = std::numeric_limits<RealT>::quiet_NaN()) const
      {
        Fixture<ScalarT> fixture(data, voltage);
        fixture.attachAllInputs(17.0);
        fixture.input(Ext::PE)   = std::isnan(pe) ? ipcmd * voltage : pe;
        fixture.input(Ext::QGEN) = std::isnan(qgen) ? iqcmd * voltage : qgen;
        if (!fixture.prepare(iqcmd, ipcmd))
        {
          return false;
        }

        auto* y  = fixture.reecb.y().getData();
        auto* yp = fixture.reecb.yp().getData();
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          y[row]  = 0.125 + 0.01 * static_cast<RealT>(row);
          yp[row] = -0.25 - 0.01 * static_cast<RealT>(row);
        }
        fixture.setCommands(iqcmd, ipcmd);
        fixture.reecb.y().setDataUpdated();
        fixture.reecb.yp().setDataUpdated();

        const auto                             y_before   = snapshot(fixture.reecb.y());
        const auto                             yp_before  = snapshot(fixture.reecb.yp());
        const auto                             bus_before = snapshot(fixture.bus.y());
        std::array<RealT, index(Ext::MAXIMUM)> input_before{};
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
          input_before[port] = fixture.input(static_cast<Ext>(port));
        }

        if (fixture.reecb.initialize() == 0)
        {
          std::cout << "Expected REECB initialization rejection: " << label << '\n';
          return false;
        }

        bool unchanged = vectorUnchanged(fixture.reecb.y(), y_before, "state")
                         && vectorUnchanged(fixture.reecb.yp(), yp_before, "derivative")
                         && vectorUnchanged(fixture.bus.y(), bus_before, "bus");
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
          unchanged &= (fixture.input(static_cast<Ext>(port)) == input_before[port]);
        }
        return unchanged;
      }

      template <typename VectorT>
      std::vector<RealT> snapshot(const VectorT& vector) const
      {
        const auto* values = vector.getData();
        return std::vector<RealT>(values, values + static_cast<size_t>(vector.getSize()));
      }

      template <typename VectorT>
      bool vectorUnchanged(const VectorT& vector, const std::vector<RealT>& expected, const char* label) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t row = 0; row < expected.size(); ++row)
        {
          if (values[row] != expected[row])
          {
            std::cout << "REECB " << label << " row " << row << " changed during rejected initialization\n";
            success = false;
          }
        }
        return success;
      }

      bool allResidualsZero(const ReecbT& reecb) const
      {
        bool        success = true;
        const auto* f       = reecb.getResidual().getData();
        const auto* yp      = reecb.yp().getData();
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          success &= scalarMatches(f[row], 0.0, "steady residual");
          success &= scalarMatches(yp[row], 0.0, "steady derivative");
        }
        return success;
      }

      bool scalarMatches(RealT actual, RealT expected, const char* label) const
      {
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "REECB " << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      void noteExpectedLogs(const char* message) const
      {
        const auto verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << message << '\n';
        Log::setVerbosity(verbosity);
      }

      template <typename T>
      void setJacobianState(Fixture<T>& fixture, RealT ilmax = 1.2) const
      {
        fixture.input(Ext::PE)     = 0.25;
        fixture.input(Ext::QGEN)   = 0.5;
        fixture.input(Ext::QEXT)   = 0.5;
        fixture.input(Ext::PFAREF) = std::atan(2.0);
        fixture.input(Ext::PREF)   = 0.25;

        auto* y               = fixture.reecb.y().getData();
        y[index(Vars::VMEAS)] = 1.0;
        y[index(Vars::PMEAS)] = 0.5;
        y[index(Vars::XPIQ)]  = 1.0;
        y[index(Vars::XPIV)]  = 0.0;
        y[index(Vars::QV)]    = 0.0;
        y[index(Vars::PORD)]  = 0.5;
        y[index(Vars::VT)]    = 1.0;
        y[index(Vars::ILMAX)] = ilmax;
        y[index(Vars::IQCMD)] = 0.1;
        y[index(Vars::IPCMD)] = 0.2;
        fixture.reecb.yp().setToConst(static_cast<T>(0.0));
        fixture.reecb.y().setDataUpdated();
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture, RealT alpha) const
      {
        auto* y     = fixture.reecb.y().getData();
        auto* yp    = fixture.reecb.yp().getData();
        auto* bus_y = fixture.bus.y().getData();
        for (size_t row = 0; row < index(Vars::MAXIMUM); ++row)
        {
          y[row].setVariableNumber(row);
          yp[row].setVariableNumber(row);
          yp[row].scaleDependencies(alpha);
        }
        for (size_t row = 0; row < static_cast<size_t>(fixture.bus.size()); ++row)
        {
          bus_y[row].setVariableNumber(index(Vars::MAXIMUM) + row);
        }
        for (size_t port = 0; port < index(Ext::MAXIMUM); ++port)
        {
          fixture.input(static_cast<Ext>(port)).setVariableNumber(fixture.inputIndex(static_cast<Ext>(port)));
        }
        fixture.reecb.y().setDataUpdated();
        fixture.reecb.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<DependencyMap> dependencyTrackingJacobian(const Data& data, RealT alpha, TestStatus& success, RealT ilmax = 1.2) const
      {
        using DepVar = DependencyTracking::Variable;
        Fixture<DepVar> fixture(data, 0.9, 0.4);
        fixture.attachAllInputs();
        fixture.input(Ext::PE)  = 0.2;
        success                *= fixture.initialize(0.0, 0.2);
        setJacobianState(fixture, ilmax);
        numberVariables(fixture, alpha);
        success *= (fixture.evaluate() == 0);

        std::vector<DependencyMap> jacobian(index(Vars::MAXIMUM));
        const auto*                residual = fixture.reecb.getResidual().getData();
        for (size_t row = 0; row < jacobian.size(); ++row)
        {
          jacobian[row] = residual[row].getDependencies();
        }
        return jacobian;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      std::vector<DependencyMap> enzymeJacobian(const Data& data, RealT alpha, TestStatus& success, RealT ilmax = 1.2) const
      {
        Fixture<ScalarT> fixture(data, 0.9, 0.4);
        fixture.attachAllInputs();
        fixture.input(Ext::PE)  = 0.2;
        success                *= fixture.initialize(0.0, 0.2);
        for (IdxT row = 0; row < fixture.bus.size(); ++row)
        {
          fixture.bus.setVariableIndex(row, fixture.reecb.size() + row);
        }
        setJacobianState(fixture, ilmax);
        fixture.reecb.updateTime(0.0, alpha);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.reecb.evaluateJacobian() == 0);
        success *= (fixture.reecb.constructCsr() == 0);
        return MapFromCsr(fixture.reecb.getCsrJacobian());
      }
#endif

      static RealT derivative(const std::vector<DependencyMap>& jacobian, size_t row, size_t column)
      {
        const auto entry = jacobian[row].find(column);
        return entry == jacobian[row].end() ? 0.0 : entry->second;
      }

      bool derivativeMatches(const std::vector<DependencyMap>& jacobian, Vars row, Vars column, RealT expected, const char* label) const
      {
        return derivativeMatches(jacobian, row, index(column), expected, label);
      }

      bool derivativeMatches(const std::vector<DependencyMap>& jacobian, Vars row, size_t column, RealT expected, const char* label) const
      {
        const RealT actual = derivative(jacobian, index(row), column);
        if (isEqual(actual, expected, kTol))
        {
          return true;
        }
        std::cout << "REECB Jacobian " << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      bool jacobiansMatch(const std::vector<DependencyMap>& dependency, const std::vector<DependencyMap>& enzyme, size_t columns) const
      {
        if (dependency.size() != enzyme.size())
        {
          std::cout << "REECB Jacobian row-count mismatch\n";
          return false;
        }

        bool success = true;
        for (size_t row = 0; row < dependency.size(); ++row)
        {
          for (size_t column = 0; column < columns; ++column)
          {
            const RealT expected = derivative(dependency, row, column);
            const RealT actual   = derivative(enzyme, row, column);
            if (!isEqual(actual, expected, kTol))
            {
              std::cout << "REECB Jacobian (" << row << ", " << column << ") backend mismatch: "
                        << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                        << actual << " != " << expected << '\n';
              success = false;
            }
          }
        }
        return success;
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
