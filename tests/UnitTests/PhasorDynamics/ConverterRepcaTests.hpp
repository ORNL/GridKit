#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ConverterRepcaTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ScalarT ir{0.0};
        ScalarT ii{0.0};
        ScalarT p{0.4};
        ScalarT q{0.1};
        ScalarT freq{1.0};
        IdxT    ir_i = 10;
        IdxT    ii_i = 11;
        IdxT    p_i  = 12;
        IdxT    q_i  = 13;
        IdxT    f_i  = 14;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> p_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> q_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freq_node;
        ir_node.set(&ir, &ir_i);
        ii_node.set(&ii, &ii_i);
        p_node.set(&p, &p_i);
        q_node.set(&q, &q_i);
        freq_node.set(&freq, &f_i);

        auto                                            data = makeData();
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        auto&                                           signals = repca.getSignals();
        signals.template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        signals.template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        signals.template attachSignalNode<Ext::PBRANCH>(&p_node);
        signals.template attachSignalNode<Ext::QBRANCH>(&q_node);
        signals.template attachSignalNode<Ext::FREQ>(&freq_node);

        success *= (repca.size() == static_cast<IdxT>(Vars::MAXIMUM));
        success *= (repca.getMonitor() != nullptr);
        success *= (repca.verify() == 0);

        auto bad                    = data;
        bad.parameters[Params::mva] = static_cast<RealT>(0.0);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> bad_mva(&bus, bad);
        bad_mva.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        bad_mva.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        bad_mva.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        bad_mva.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        bad_mva.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        success *= (bad_mva.verify() > 0);

        bad                         = data;
        bad.parameters[Params::Tfv] = static_cast<RealT>(0.0);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> bad_tfv(&bus, bad);
        bad_tfv.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        bad_tfv.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        bad_tfv.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        bad_tfv.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        bad_tfv.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        success *= (bad_tfv.verify() > 0);

        bad                          = data;
        bad.parameters[Params::Qmin] = static_cast<RealT>(2.0);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> bad_q_limits(&bus, bad);
        bad_q_limits.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        bad_q_limits.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        bad_q_limits.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        bad_q_limits.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        bad_q_limits.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        success *= (bad_q_limits.verify() > 0);

        bad                           = data;
        bad.parameters[Params::fdbd1] = static_cast<RealT>(0.1);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> bad_deadband(&bus, bad);
        bad_deadband.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        bad_deadband.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        bad_deadband.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        bad_deadband.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        bad_deadband.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        success *= (bad_deadband.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome signalVerification()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        success *= (repca.verify() > 0);

        ScalarT ir{0.0};
        ScalarT ii{0.0};
        ScalarT p{0.4};
        ScalarT q{0.1};
        ScalarT freq{1.0};
        ScalarT vref{1.0};
        IdxT    ir_i = 20;
        IdxT    ii_i = 21;
        IdxT    p_i  = 22;
        IdxT    q_i  = 23;
        IdxT    f_i  = 24;
        IdxT    v_i  = 25;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> p_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> q_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freq_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vref_node;
        ir_node.set(&ir, &ir_i);
        ii_node.set(&ii, &ii_i);
        p_node.set(&p, &p_i);
        q_node.set(&q, &q_i);
        freq_node.set(&freq, &f_i);

        repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        repca.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        success *= (repca.verify() == 0);

        repca.getSignals().template attachSignalNode<Ext::VREF>(&vref_node);
        success *= (repca.verify() > 0);
        vref_node.set(&vref, &v_i);
        success *= (repca.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome initializationAndResidual()
      {
        TestStatus success = true;

        auto data                         = makeData();
        data.parameters[Params::mva]      = static_cast<RealT>(50.0);
        data.parameters[Params::Freqflag] = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        ScalarT ir{0.2};
        ScalarT ii{-0.1};
        ScalarT p{0.4};
        ScalarT q{0.1};
        ScalarT freq{1.0};
        ScalarT freqref{-4.0};
        ScalarT vref{-3.0};
        ScalarT qref{-2.0};
        ScalarT pplantref{-1.0};
        IdxT    ir_i       = 30;
        IdxT    ii_i       = 31;
        IdxT    p_i        = 32;
        IdxT    q_i        = 33;
        IdxT    f_i        = 34;
        IdxT    freqref_i  = 35;
        IdxT    vref_i     = 36;
        IdxT    qref_i     = 37;
        IdxT    plantref_i = 38;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> p_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> q_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freq_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freqref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> plantref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qext_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pext_node;

        ir_node.set(&ir, &ir_i);
        ii_node.set(&ii, &ii_i);
        p_node.set(&p, &p_i);
        q_node.set(&q, &q_i);
        freq_node.set(&freq, &f_i);
        freqref_node.set(&freqref, &freqref_i);
        vref_node.set(&vref, &vref_i);
        qref_node.set(&qref, &qref_i);
        plantref_node.set(&pplantref, &plantref_i);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        auto&                                           signals = repca.getSignals();
        signals.template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        signals.template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        signals.template attachSignalNode<Ext::PBRANCH>(&p_node);
        signals.template attachSignalNode<Ext::QBRANCH>(&q_node);
        signals.template attachSignalNode<Ext::FREQ>(&freq_node);
        signals.template attachSignalNode<Ext::FREQREF>(&freqref_node);
        signals.template attachSignalNode<Ext::VREF>(&vref_node);
        signals.template attachSignalNode<Ext::QREF>(&qref_node);
        signals.template attachSignalNode<Ext::PPLANTREF>(&plantref_node);
        signals.template assignSignalNode<Vars::QEXT>(&qext_node);
        signals.template assignSignalNode<Vars::PEXT>(&pext_node);

        success                                *= (repca.allocate() == 0);
        repca.y().getData()[index(Vars::QEXT)]  = static_cast<ScalarT>(0.25);
        repca.y().getData()[index(Vars::PEXT)]  = static_cast<ScalarT>(0.70);
        repca.y().setDataUpdated();
        success *= (repca.initialize() == 0);
        success *= (repca.tagDifferentiable() == 0);
        success *= (repca.evaluateResidual() == 0);

        const ScalarT base_ratio = static_cast<ScalarT>(2.0);
        const ScalarT qmeas      = base_ratio * q;
        const ScalarT pmeas      = base_ratio * p;
        const ScalarT pfreq0 =
            static_cast<RealT>(20.0) * Math::ramp(static_cast<ScalarT>(0.0));

        success *= isEqual(repca.y().getData()[index(Vars::QEXT)], static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(repca.y().getData()[index(Vars::PEXT)], static_cast<ScalarT>(0.70), kTol);
        success *= isEqual(qext_node.read(), static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(pext_node.read(), static_cast<ScalarT>(0.70), kTol);
        success *= isEqual(freqref, freq, kTol);
        success *= isEqual(vref, repca.y().getData()[index(Vars::VMEAS)], kTol);
        success *= isEqual(qref, qmeas / base_ratio, kTol);
        success *= isEqual(pplantref, (pmeas - pfreq0) / base_ratio, kTol);

        for (size_t i = 0; i < repca.getResidual().getSize(); ++i)
        {
          if (!isEqual(repca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "REPCA residual row " << i << " is "
                      << repca.getResidual().getData()[i] << "\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

      TestOutcome residualEquations()
      {
        TestStatus success = true;

        auto                               data = makeResidualData();
        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.9, 0.4);
        bus.allocate();
        bus.initialize();

        ScalarT ir{0.08};
        ScalarT ii{-0.02};
        ScalarT p{0.41};
        ScalarT q{0.13};
        ScalarT freq{0.99};
        ScalarT freqref{1.0};
        ScalarT vref{1.01};
        ScalarT qref{0.12};
        ScalarT pplantref{0.55};
        IdxT    ir_i       = 40;
        IdxT    ii_i       = 41;
        IdxT    p_i        = 42;
        IdxT    q_i        = 43;
        IdxT    f_i        = 44;
        IdxT    freqref_i  = 45;
        IdxT    vref_i     = 46;
        IdxT    qref_i     = 47;
        IdxT    plantref_i = 48;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> p_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> q_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freq_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freqref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> plantref_node;

        ir_node.set(&ir, &ir_i);
        ii_node.set(&ii, &ii_i);
        p_node.set(&p, &p_i);
        q_node.set(&q, &q_i);
        freq_node.set(&freq, &f_i);
        freqref_node.set(&freqref, &freqref_i);
        vref_node.set(&vref, &vref_i);
        qref_node.set(&qref, &qref_i);
        plantref_node.set(&pplantref, &plantref_i);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        repca.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        repca.getSignals().template attachSignalNode<Ext::FREQREF>(&freqref_node);
        repca.getSignals().template attachSignalNode<Ext::VREF>(&vref_node);
        repca.getSignals().template attachSignalNode<Ext::QREF>(&qref_node);
        repca.getSignals().template attachSignalNode<Ext::PPLANTREF>(&plantref_node);

        success *= (repca.allocate() == 0);

        repca.y().getData()[index(Vars::VMEAS)]  = static_cast<ScalarT>(0.98);
        repca.y().getData()[index(Vars::QMEAS)]  = static_cast<ScalarT>(0.11);
        repca.y().getData()[index(Vars::XQPI)]   = static_cast<ScalarT>(0.07);
        repca.y().getData()[index(Vars::XQLAG)]  = static_cast<ScalarT>(0.14);
        repca.y().getData()[index(Vars::PMEAS)]  = static_cast<ScalarT>(0.44);
        repca.y().getData()[index(Vars::XPPI)]   = static_cast<ScalarT>(0.21);
        repca.y().getData()[index(Vars::PREF)]   = static_cast<ScalarT>(0.60);
        repca.y().getData()[index(Vars::V)]      = static_cast<ScalarT>(1.0);
        repca.y().getData()[index(Vars::VLDC)]   = static_cast<ScalarT>(0.99);
        repca.y().getData()[index(Vars::VDROOP)] = static_cast<ScalarT>(1.05);
        repca.y().getData()[index(Vars::VCTRL)]  = static_cast<ScalarT>(1.02);
        repca.y().getData()[index(Vars::SFRZ)]   = static_cast<ScalarT>(0.8);
        repca.y().getData()[index(Vars::ERQ)]    = static_cast<ScalarT>(0.03);
        repca.y().getData()[index(Vars::ERQDB)]  = static_cast<ScalarT>(0.02);
        repca.y().getData()[index(Vars::ERQLIM)] = static_cast<ScalarT>(0.02);
        repca.y().getData()[index(Vars::QPI)]    = static_cast<ScalarT>(0.27);
        repca.y().getData()[index(Vars::QEXT)]   = static_cast<ScalarT>(0.20);
        repca.y().getData()[index(Vars::EF)]     = static_cast<ScalarT>(0.01);
        repca.y().getData()[index(Vars::EP)]     = static_cast<ScalarT>(0.04);
        repca.y().getData()[index(Vars::EPLIM)]  = static_cast<ScalarT>(0.04);
        repca.y().getData()[index(Vars::PPI)]    = static_cast<ScalarT>(0.66);
        repca.y().getData()[index(Vars::PEXT)]   = static_cast<ScalarT>(0.61);

        repca.yp().getData()[index(Vars::VMEAS)] = static_cast<ScalarT>(0.01);
        repca.yp().getData()[index(Vars::QMEAS)] = static_cast<ScalarT>(-0.02);
        repca.yp().getData()[index(Vars::XQPI)]  = static_cast<ScalarT>(0.03);
        repca.yp().getData()[index(Vars::XQLAG)] = static_cast<ScalarT>(-0.04);
        repca.yp().getData()[index(Vars::PMEAS)] = static_cast<ScalarT>(0.02);
        repca.yp().getData()[index(Vars::XPPI)]  = static_cast<ScalarT>(-0.01);
        repca.yp().getData()[index(Vars::PREF)]  = static_cast<ScalarT>(0.05);

        repca.y().setDataUpdated();
        repca.yp().setDataUpdated();

        success *= (repca.evaluateResidual() == 0);

        const ScalarT vldc_r = bus.Vr() - static_cast<RealT>(0.02) * ir
                               + static_cast<RealT>(0.03) * ii;
        const ScalarT vldc_i = bus.Vi() - static_cast<RealT>(0.02) * ii
                               - static_cast<RealT>(0.03) * ir;
        const ScalarT pfreq = static_cast<RealT>(2.0) * Math::ramp(repca.y().getData()[index(Vars::EF)])
                              - Math::ramp(-repca.y().getData()[index(Vars::EF)]);

        std::vector<ScalarT> expected(static_cast<size_t>(Vars::MAXIMUM),
                                      static_cast<ScalarT>(0.0));
        expected[index(Vars::VMEAS)] =
            -repca.yp().getData()[index(Vars::VMEAS)]
            + (repca.y().getData()[index(Vars::VCTRL)] - repca.y().getData()[index(Vars::VMEAS)])
                  / static_cast<RealT>(0.2);
        expected[index(Vars::QMEAS)] =
            -repca.yp().getData()[index(Vars::QMEAS)]
            + (q - repca.y().getData()[index(Vars::QMEAS)]) / static_cast<RealT>(0.2);
        expected[index(Vars::XQPI)] =
            -repca.yp().getData()[index(Vars::XQPI)]
            + repca.y().getData()[index(Vars::SFRZ)]
                  * Math::antiwindup(repca.y().getData()[index(Vars::QPI)],
                                     static_cast<RealT>(3.0)
                                         * repca.y().getData()[index(Vars::ERQLIM)],
                                     static_cast<RealT>(-0.8),
                                     static_cast<RealT>(0.9));
        expected[index(Vars::XQLAG)] =
            -repca.yp().getData()[index(Vars::XQLAG)]
            + (repca.y().getData()[index(Vars::QPI)] - repca.y().getData()[index(Vars::XQLAG)])
                  / static_cast<RealT>(1.5);
        expected[index(Vars::PMEAS)] =
            -repca.yp().getData()[index(Vars::PMEAS)]
            + (p - repca.y().getData()[index(Vars::PMEAS)]) / static_cast<RealT>(0.4);
        expected[index(Vars::XPPI)] =
            -repca.yp().getData()[index(Vars::XPPI)]
            + Math::antiwindup(repca.y().getData()[index(Vars::PPI)],
                               static_cast<RealT>(1.8) * repca.y().getData()[index(Vars::EPLIM)],
                               static_cast<RealT>(0.0),
                               static_cast<RealT>(1.2));
        expected[index(Vars::PREF)] =
            -repca.yp().getData()[index(Vars::PREF)]
            + (repca.y().getData()[index(Vars::PPI)] - repca.y().getData()[index(Vars::PREF)])
                  / static_cast<RealT>(0.5);

        expected[index(Vars::V)] =
            -repca.y().getData()[index(Vars::V)] * repca.y().getData()[index(Vars::V)]
            + bus.Vr() * bus.Vr() + bus.Vi() * bus.Vi();
        expected[index(Vars::VLDC)] =
            -repca.y().getData()[index(Vars::VLDC)] * repca.y().getData()[index(Vars::VLDC)]
            + vldc_r * vldc_r + vldc_i * vldc_i;
        expected[index(Vars::VDROOP)] =
            -repca.y().getData()[index(Vars::VDROOP)] + repca.y().getData()[index(Vars::V)]
            + static_cast<RealT>(0.4) * q;
        expected[index(Vars::VCTRL)] =
            -repca.y().getData()[index(Vars::VCTRL)] + repca.y().getData()[index(Vars::VLDC)];
        expected[index(Vars::SFRZ)] =
            -repca.y().getData()[index(Vars::SFRZ)]
            + Math::above(repca.y().getData()[index(Vars::V)], static_cast<RealT>(0.7));
        expected[index(Vars::ERQ)] =
            -repca.y().getData()[index(Vars::ERQ)] + vref - repca.y().getData()[index(Vars::VMEAS)];
        expected[index(Vars::ERQDB)] =
            -repca.y().getData()[index(Vars::ERQDB)]
            + Math::deadband2(repca.y().getData()[index(Vars::ERQ)],
                              static_cast<RealT>(-0.02),
                              static_cast<RealT>(0.03));
        expected[index(Vars::ERQLIM)] =
            -repca.y().getData()[index(Vars::ERQLIM)]
            + Math::clamp(repca.y().getData()[index(Vars::ERQDB)],
                          static_cast<RealT>(-0.7),
                          static_cast<RealT>(0.8));
        expected[index(Vars::QPI)] =
            -repca.y().getData()[index(Vars::QPI)]
            + Math::clamp(static_cast<RealT>(2.0) * repca.y().getData()[index(Vars::ERQLIM)]
                              + repca.y().getData()[index(Vars::XQPI)],
                          static_cast<RealT>(-0.8),
                          static_cast<RealT>(0.9));
        expected[index(Vars::QEXT)] =
            -static_cast<RealT>(1.5)
                * (repca.y().getData()[index(Vars::QEXT)] - repca.y().getData()[index(Vars::XQLAG)])
            + static_cast<RealT>(0.2)
                  * (repca.y().getData()[index(Vars::QPI)] - repca.y().getData()[index(Vars::XQLAG)]);

        expected[index(Vars::EF)] =
            -repca.y().getData()[index(Vars::EF)]
            + Math::deadband2(freqref - freq,
                              static_cast<RealT>(-0.01),
                              static_cast<RealT>(0.015));
        expected[index(Vars::EP)] =
            -repca.y().getData()[index(Vars::EP)] + pplantref - repca.y().getData()[index(Vars::PMEAS)]
            + pfreq;
        expected[index(Vars::EPLIM)] =
            -repca.y().getData()[index(Vars::EPLIM)]
            + Math::clamp(repca.y().getData()[index(Vars::EP)],
                          static_cast<RealT>(-0.5),
                          static_cast<RealT>(0.6));
        expected[index(Vars::PPI)] =
            -repca.y().getData()[index(Vars::PPI)]
            + Math::clamp(static_cast<RealT>(1.7) * repca.y().getData()[index(Vars::EPLIM)]
                              + repca.y().getData()[index(Vars::XPPI)],
                          static_cast<RealT>(0.0),
                          static_cast<RealT>(1.2));
        expected[index(Vars::PEXT)] =
            -repca.y().getData()[index(Vars::PEXT)] + repca.y().getData()[index(Vars::PREF)];

        for (size_t i = 0; i < expected.size(); ++i)
        {
          if (!isEqual(repca.getResidual().getData()[i], expected[i], kTol))
          {
            std::cout << "REPCA residual mismatch at row " << i << ": "
                      << repca.getResidual().getData()[i] << " != " << expected[i] << "\n";
            success = false;
          }
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
    "case_name": "REPCA parser",
    "case_description": "REPCA parser behavior test",
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
  "signals": [
    { "signal_id": 10, "name": "Ir" },
    { "signal_id": 11, "name": "Ii" },
    { "signal_id": 12, "name": "Pbranch" },
    { "signal_id": 13, "name": "Qbranch" },
    { "signal_id": 14, "name": "Freq" },
    { "signal_id": 15, "name": "Qext" },
    { "signal_id": 16, "name": "Pext" }
  ],
  "devices": [
    {
      "class": "Repca",
      "ports": {
        "bus": 1,
        "ibranchr": 10,
        "ibranchi": 11,
        "pbranch": 12,
        "qbranch": 13,
        "freq": 14,
        "qext": 15,
        "pext": 16
      },
      "id": "RP1",
      "params": { "mva": 100.0 },
      "mon": ["qext", "pext", "vmeas", "qmeas", "pmeas"]
    }
  ]
}
)json");

        auto data  = PhasorDynamics::parseSystemModelData(input);
        success   *= (data.repca.size() == 1);
        success   *= (data.repca[0].device_class == "Repca");
        success   *= (data.repca[0].buses.at(PhasorDynamics::Converter::RepcaBuses::bus) == 1);
        success   *= (data.repca[0].signal_inputs.at(PhasorDynamics::Converter::RepcaSignalInputs::ibranchr) == 10);
        success   *= (data.repca[0].signal_inputs.at(PhasorDynamics::Converter::RepcaSignalInputs::ibranchi) == 11);
        success   *= (data.repca[0].signal_inputs.at(PhasorDynamics::Converter::RepcaSignalInputs::pbranch) == 12);
        success   *= (data.repca[0].signal_inputs.at(PhasorDynamics::Converter::RepcaSignalInputs::qbranch) == 13);
        success   *= (data.repca[0].signal_inputs.at(PhasorDynamics::Converter::RepcaSignalInputs::freq) == 14);
        success   *= (data.repca[0].signal_outputs.at(PhasorDynamics::Converter::RepcaSignalOutputs::qext) == 15);
        success   *= (data.repca[0].signal_outputs.at(PhasorDynamics::Converter::RepcaSignalOutputs::pext) == 16);
        success   *= (std::get_if<double>(
                        &data.repca[0].parameters.at(Params::mva))
                    != nullptr);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.size() == 0);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto dependency_tracking_jacobian = dependencyTrackingJacobian();
        auto enzyme_jacobian              = enzymeJacobian();

        success *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
        const auto nrows =
            std::min(dependency_tracking_jacobian.size(), enzyme_jacobian.size());
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
      using Params = PhasorDynamics::Converter::RepcaParameters;
      using Vars   = PhasorDynamics::Converter::RepcaInternalVariables;
      using Ext    = PhasorDynamics::Converter::RepcaExternalVariables;

      static size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      auto makeData() -> PhasorDynamics::Converter::RepcaData<RealT, IdxT>
      {
        using Mon = PhasorDynamics::Converter::RepcaMonitorableVariables;

        PhasorDynamics::Converter::RepcaData<RealT, IdxT> data;
        data.device_class          = "Repca";
        data.disambiguation_string = "repca_test";
        data.monitored_variables.insert(Mon::qext);
        data.monitored_variables.insert(Mon::pext);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::qmeas);
        data.monitored_variables.insert(Mon::pmeas);
        return data;
      }

      auto makeResidualData() -> PhasorDynamics::Converter::RepcaData<RealT, IdxT>
      {
        auto data = makeData();

        data.parameters[Params::mva]       = static_cast<RealT>(100.0);
        data.parameters[Params::VcompFlag] = true;
        data.parameters[Params::RefFlag]   = true;
        data.parameters[Params::Freqflag]  = true;
        data.parameters[Params::Tfltr]     = static_cast<RealT>(0.2);
        data.parameters[Params::Vfrz]      = static_cast<RealT>(0.7);
        data.parameters[Params::Rc]        = static_cast<RealT>(0.02);
        data.parameters[Params::Xc]        = static_cast<RealT>(0.03);
        data.parameters[Params::Kc]        = static_cast<RealT>(0.4);
        data.parameters[Params::dbdlow]    = static_cast<RealT>(-0.02);
        data.parameters[Params::dbdupper]  = static_cast<RealT>(0.03);
        data.parameters[Params::emax]      = static_cast<RealT>(0.8);
        data.parameters[Params::emin]      = static_cast<RealT>(-0.7);
        data.parameters[Params::Kp]        = static_cast<RealT>(2.0);
        data.parameters[Params::Ki]        = static_cast<RealT>(3.0);
        data.parameters[Params::Qmax]      = static_cast<RealT>(0.9);
        data.parameters[Params::Qmin]      = static_cast<RealT>(-0.8);
        data.parameters[Params::Tft]       = static_cast<RealT>(0.2);
        data.parameters[Params::Tfv]       = static_cast<RealT>(1.5);
        data.parameters[Params::Tp]        = static_cast<RealT>(0.4);
        data.parameters[Params::fdbd1]     = static_cast<RealT>(-0.01);
        data.parameters[Params::fdbd2]     = static_cast<RealT>(0.015);
        data.parameters[Params::Ddn]       = static_cast<RealT>(2.0);
        data.parameters[Params::Dup]       = static_cast<RealT>(1.0);
        data.parameters[Params::femax]     = static_cast<RealT>(0.6);
        data.parameters[Params::femin]     = static_cast<RealT>(-0.5);
        data.parameters[Params::Kpg]       = static_cast<RealT>(1.7);
        data.parameters[Params::Kig]       = static_cast<RealT>(1.8);
        data.parameters[Params::Pmax]      = static_cast<RealT>(1.2);
        data.parameters[Params::Pmin]      = static_cast<RealT>(0.0);
        data.parameters[Params::Tlag]      = static_cast<RealT>(0.5);
        return data;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      void setValue(double& target, double value)
      {
        target = value;
      }

      void setValue(DependencyTracking::Variable& target, double value)
      {
        target.setValue(value);
      }

      template <typename model_type, typename bus_type>
      void setJacobianState(model_type& repca, bus_type& bus)
      {
        setValue(bus.Vr(), 0.9);
        setValue(bus.Vi(), 0.4);

        setValue(repca.y().getData()[index(Vars::VMEAS)], 0.98);
        setValue(repca.y().getData()[index(Vars::QMEAS)], 0.11);
        setValue(repca.y().getData()[index(Vars::XQPI)], 0.07);
        setValue(repca.y().getData()[index(Vars::XQLAG)], 0.14);
        setValue(repca.y().getData()[index(Vars::PMEAS)], 0.44);
        setValue(repca.y().getData()[index(Vars::XPPI)], 0.21);
        setValue(repca.y().getData()[index(Vars::PREF)], 0.60);
        setValue(repca.y().getData()[index(Vars::V)], 1.0);
        setValue(repca.y().getData()[index(Vars::VLDC)], 0.99);
        setValue(repca.y().getData()[index(Vars::VDROOP)], 1.05);
        setValue(repca.y().getData()[index(Vars::VCTRL)], 1.02);
        setValue(repca.y().getData()[index(Vars::SFRZ)], 0.8);
        setValue(repca.y().getData()[index(Vars::ERQ)], 0.03);
        setValue(repca.y().getData()[index(Vars::ERQDB)], 0.02);
        setValue(repca.y().getData()[index(Vars::ERQLIM)], 0.02);
        setValue(repca.y().getData()[index(Vars::QPI)], 0.27);
        setValue(repca.y().getData()[index(Vars::QEXT)], 0.20);
        setValue(repca.y().getData()[index(Vars::EF)], 0.01);
        setValue(repca.y().getData()[index(Vars::EP)], 0.04);
        setValue(repca.y().getData()[index(Vars::EPLIM)], 0.04);
        setValue(repca.y().getData()[index(Vars::PPI)], 0.66);
        setValue(repca.y().getData()[index(Vars::PEXT)], 0.61);

        setValue(repca.yp().getData()[index(Vars::VMEAS)], 0.01);
        setValue(repca.yp().getData()[index(Vars::QMEAS)], -0.02);
        setValue(repca.yp().getData()[index(Vars::XQPI)], 0.03);
        setValue(repca.yp().getData()[index(Vars::XQLAG)], -0.04);
        setValue(repca.yp().getData()[index(Vars::PMEAS)], 0.02);
        setValue(repca.yp().getData()[index(Vars::XPPI)], -0.01);
        setValue(repca.yp().getData()[index(Vars::PREF)], 0.05);

        bus.y().setDataUpdated();
        repca.y().setDataUpdated();
        repca.yp().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        auto data = makeResidualData();

        PhasorDynamics::Bus<DepVar, IdxT>              bus(DepVar{0.9}, DepVar{0.4});
        PhasorDynamics::Converter::Repca<DepVar, IdxT> repca(&bus, data);

        DepVar ir{0.08};
        DepVar ii{-0.02};
        DepVar p{0.41};
        DepVar q{0.13};
        DepVar freq{0.99};
        DepVar freqref{1.0};
        DepVar vref{1.01};
        DepVar qref{0.12};
        DepVar pplantref{0.55};
        IdxT   ir_i       = 24;
        IdxT   ii_i       = 25;
        IdxT   p_i        = 26;
        IdxT   q_i        = 27;
        IdxT   f_i        = 28;
        IdxT   freqref_i  = 29;
        IdxT   vref_i     = 30;
        IdxT   qref_i     = 31;
        IdxT   plantref_i = 32;

        ir.setVariableNumber(ir_i);
        ii.setVariableNumber(ii_i);
        p.setVariableNumber(p_i);
        q.setVariableNumber(q_i);
        freq.setVariableNumber(f_i);
        freqref.setVariableNumber(freqref_i);
        vref.setVariableNumber(vref_i);
        qref.setVariableNumber(qref_i);
        pplantref.setVariableNumber(plantref_i);

        PhasorDynamics::SignalNode<DepVar, IdxT> ir_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> ii_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> p_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> q_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> freq_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> freqref_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> vref_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> qref_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> plantref_node;
        ir_node.set(&ir, &ir_i);
        ii_node.set(&ii, &ii_i);
        p_node.set(&p, &p_i);
        q_node.set(&q, &q_i);
        freq_node.set(&freq, &f_i);
        freqref_node.set(&freqref, &freqref_i);
        vref_node.set(&vref, &vref_i);
        qref_node.set(&qref, &qref_i);
        plantref_node.set(&pplantref, &plantref_i);

        repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        repca.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        repca.getSignals().template attachSignalNode<Ext::FREQREF>(&freqref_node);
        repca.getSignals().template attachSignalNode<Ext::VREF>(&vref_node);
        repca.getSignals().template attachSignalNode<Ext::QREF>(&qref_node);
        repca.getSignals().template attachSignalNode<Ext::PPLANTREF>(&plantref_node);

        bus.allocate();
        repca.allocate();
        for (IdxT i = 0; i < repca.size(); ++i)
        {
          repca.y().getData()[static_cast<size_t>(i)].setVariableNumber(i);
          repca.yp().getData()[static_cast<size_t>(i)].setVariableNumber(i);
        }
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.y().getData()[static_cast<size_t>(i)].setVariableNumber(i + repca.size());
        }

        repca.y().setDataUpdated();
        repca.yp().setDataUpdated();
        bus.y().setDataUpdated();

        setJacobianState(repca, bus);
        repca.evaluateResidual();

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(
            static_cast<size_t>(repca.size()));
        for (IdxT i = 0; i < repca.size(); ++i)
        {
          dependencies[static_cast<size_t>(i)] =
              repca.getResidual().getData()[static_cast<size_t>(i)].getDependencies();
        }
        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian()
      {
        auto data = makeResidualData();

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(0.9, 0.4);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);

        ScalarT ir{0.08};
        ScalarT ii{-0.02};
        ScalarT p{0.41};
        ScalarT q{0.13};
        ScalarT freq{0.99};
        ScalarT freqref{1.0};
        ScalarT vref{1.01};
        ScalarT qref{0.12};
        ScalarT pplantref{0.55};
        IdxT    ir_i       = 24;
        IdxT    ii_i       = 25;
        IdxT    p_i        = 26;
        IdxT    q_i        = 27;
        IdxT    f_i        = 28;
        IdxT    freqref_i  = 29;
        IdxT    vref_i     = 30;
        IdxT    qref_i     = 31;
        IdxT    plantref_i = 32;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ir_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ii_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> p_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> q_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freq_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> freqref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> plantref_node;
        ir_node.set(&ir, &ir_i);
        ii_node.set(&ii, &ii_i);
        p_node.set(&p, &p_i);
        q_node.set(&q, &q_i);
        freq_node.set(&freq, &f_i);
        freqref_node.set(&freqref, &freqref_i);
        vref_node.set(&vref, &vref_i);
        qref_node.set(&qref, &qref_i);
        plantref_node.set(&pplantref, &plantref_i);

        repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ir_node);
        repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ii_node);
        repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&p_node);
        repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&q_node);
        repca.getSignals().template attachSignalNode<Ext::FREQ>(&freq_node);
        repca.getSignals().template attachSignalNode<Ext::FREQREF>(&freqref_node);
        repca.getSignals().template attachSignalNode<Ext::VREF>(&vref_node);
        repca.getSignals().template attachSignalNode<Ext::QREF>(&qref_node);
        repca.getSignals().template attachSignalNode<Ext::PPLANTREF>(&plantref_node);

        bus.allocate();
        repca.allocate();
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + repca.size());
          bus.setResidualIndex(i, i + repca.size());
        }

        setJacobianState(repca, bus);
        repca.updateTime(0.0, 1.0);
        repca.evaluateResidual();
        repca.evaluateJacobian();
        repca.constructCsr();
        return MapFromCsr(repca.getCsrJacobian());
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
