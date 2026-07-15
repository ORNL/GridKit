#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ConverterReecbTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using ReecbT  = PhasorDynamics::Converter::Reecb<ScalarT, IdxT>;
      using Var     = PhasorDynamics::Converter::ReecbInternalVariables;
      using Ext     = PhasorDynamics::Converter::ReecbExternalVariables;
      using Params  = PhasorDynamics::Converter::ReecbParameters;
      using Buses   = PhasorDynamics::Converter::ReecbBuses;
      using Outputs = PhasorDynamics::Converter::ReecbSignalOutputs;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        ReecbT reecb(&bus, makeData());
        success *= (reecb.size() == static_cast<IdxT>(Var::MAXIMUM));
        success *= (reecb.getMonitor() != nullptr);
        success *= (reecb.verify() == 0);

        auto minimal                    = makeMinimalData();
        minimal.parameters[Params::mva] = static_cast<RealT>(100.0);
        ReecbT minimal_model(&bus, minimal);
        success *= (minimal_model.verify() == 0);

        ReecbT missing_mva(&bus, makeMinimalData());
        success *= (missing_mva.verify() > 0);

        auto bad_band                     = makeData();
        bad_band.parameters[Params::Vdip] = static_cast<RealT>(1.2);
        bad_band.parameters[Params::Vup]  = static_cast<RealT>(1.2);
        ReecbT bad_band_model(&bus, bad_band);
        success *= (bad_band_model.verify() > 0);

        auto bad_imax                     = makeData();
        bad_imax.parameters[Params::Imax] = static_cast<RealT>(-1.0);
        ReecbT bad_imax_model(&bus, bad_imax);
        success *= (bad_imax_model.verify() > 0);

        ScalarT                                   pe_value{0.5};
        IdxT                                      pe_index = 20;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        pe_node.set(&pe_value, &pe_index);

        ReecbT half_connected(&bus, makeData());
        half_connected.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        success *= (half_connected.verify() == 0);

        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_pe_node;
        ReecbT                                    unlinked(&bus, makeData());
        unlinked.getSignals().template attachSignalNode<Ext::PE>(&unlinked_pe_node);
        success *= (unlinked.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome signals()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.6};
        IdxT    iqcmd_index = 21;
        IdxT    ipcmd_index = 22;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        ReecbT reecb(&bus, makeData());
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.6));
        success *= (reecb.verify() == 0);
        success *= (reecb.initialize() == 0);
        success *= (reecb.tagDifferentiable() == 0);
        success *= (reecb.evaluateResidual() == 0);

        success *= isEqual(reecb.y().getData()[index(Var::VMEAS)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::PMEAS)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::QREF)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::PORD)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(iqcmd_node.read(), reecb.y().getData()[index(Var::IQCMD)], kTol);
        success *= isEqual(ipcmd_node.read(), reecb.y().getData()[index(Var::IPCMD)], kTol);
        success *= (reecb.tag()[index(Var::VMEAS)] == true);
        success *= (reecb.tag()[index(Var::PMEAS)] == true);
        success *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome publishRefs()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        auto data                     = makeData();
        data.parameters[Params::Qmin] = static_cast<RealT>(-2.0);
        data.parameters[Params::Qmax] = static_cast<RealT>(2.0);
        data.parameters[Params::Pmax] = static_cast<RealT>(2.0);
        data.parameters[Params::Imax] = static_cast<RealT>(2.0);

        ScalarT iqcmd_value{-0.2};
        ScalarT ipcmd_value{1.6};
        ScalarT qext_value{99.0};
        ScalarT pfaref_value{99.0};
        ScalarT pref_value{99.0};
        IdxT    iqcmd_index = 30;
        IdxT    ipcmd_index = 31;
        IdxT    qext_index  = 32;
        IdxT    pfref_index = 33;
        IdxT    pref_index  = 34;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qext_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pfaref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        qext_node.set(&qext_value, &qext_index);
        pfaref_node.set(&pfaref_value, &pfref_index);
        pref_node.set(&pref_value, &pref_index);

        ReecbT reecb(&bus, data);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);
        reecb.getSignals().template attachSignalNode<Ext::QEXT>(&qext_node);
        reecb.getSignals().template attachSignalNode<Ext::PFAREF>(&pfaref_node);
        reecb.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);

        success *= (reecb.allocate() == 0);
        iqcmd_node.init(static_cast<ScalarT>(-0.2));
        ipcmd_node.init(static_cast<ScalarT>(1.6));
        success *= (reecb.verify() == 0);
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        const ScalarT expected_pfaref  = static_cast<ScalarT>(std::atan(static_cast<RealT>(-0.2 / 1.6)));
        success                       *= isEqual(qext_node.read(), static_cast<ScalarT>(-0.2), kTol);
        success                       *= isEqual(pfaref_node.read(), expected_pfaref, kTol);
        success                       *= isEqual(pref_node.read(), static_cast<ScalarT>(1.6), kTol);
        success                       *= isEqual(reecb.y().getData()[index(Var::QREF)], static_cast<ScalarT>(-0.2), kTol);
        success                       *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome baseSignals()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        ScalarT pe_value{99.0};
        ScalarT qgen_value{99.0};
        ScalarT qext_value{99.0};
        ScalarT pfaref_value{99.0};
        ScalarT pref_value{99.0};
        ScalarT iqcmd_value{0.0};
        ScalarT ipcmd_value{0.0};
        IdxT    pe_index     = 40;
        IdxT    qgen_index   = 41;
        IdxT    qext_index   = 42;
        IdxT    pfaref_index = 43;
        IdxT    pref_index   = 44;
        IdxT    iqcmd_index  = 45;
        IdxT    ipcmd_index  = 46;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qext_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pfaref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        qext_node.set(&qext_value, &qext_index);
        pfaref_node.set(&pfaref_value, &pfaref_index);
        pref_node.set(&pref_value, &pref_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        ReecbT reecb(&bus, data);
        reecb.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reecb.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reecb.getSignals().template attachSignalNode<Ext::QEXT>(&qext_node);
        reecb.getSignals().template attachSignalNode<Ext::PFAREF>(&pfaref_node);
        reecb.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.05));
        ipcmd_node.init(static_cast<ScalarT>(0.25));
        success *= (reecb.verify() == 0);
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        const ScalarT expected_pfaref  = static_cast<ScalarT>(std::atan(static_cast<RealT>(0.1 / 0.5)));
        success                       *= isEqual(reecb.y().getData()[index(Var::PMEAS)], static_cast<ScalarT>(0.5), kTol);
        success                       *= isEqual(reecb.y().getData()[index(Var::QREF)], static_cast<ScalarT>(0.1), kTol);
        success                       *= isEqual(reecb.y().getData()[index(Var::PORD)], static_cast<ScalarT>(0.5), kTol);
        success                       *= isEqual(pe_node.read(), static_cast<ScalarT>(0.25), kTol);
        success                       *= isEqual(qgen_node.read(), static_cast<ScalarT>(0.05), kTol);
        success                       *= isEqual(qext_node.read(), static_cast<ScalarT>(0.05), kTol);
        success                       *= isEqual(pfaref_node.read(), expected_pfaref, kTol);
        success                       *= isEqual(pref_node.read(), static_cast<ScalarT>(0.25), kTol);
        success                       *= isEqual(iqcmd_node.read(), static_cast<ScalarT>(0.05), kTol);
        success                       *= isEqual(ipcmd_node.read(), static_cast<ScalarT>(0.25), kTol);
        success                       *= isEqual(reecb.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.05), kTol);
        success                       *= isEqual(reecb.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.25), kTol);
        success                       *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome feedbackBase()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        ScalarT pe_value{0.25};
        ScalarT qgen_value{0.05};
        ScalarT iqcmd_value{0.0};
        ScalarT ipcmd_value{0.0};
        IdxT    pe_index    = 40;
        IdxT    qgen_index  = 41;
        IdxT    iqcmd_index = 42;
        IdxT    ipcmd_index = 43;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        ReecbT reecb(&bus, data);
        reecb.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reecb.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.05));
        ipcmd_node.init(static_cast<ScalarT>(0.25));
        success *= (reecb.verify() == 0);
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        success *= isEqual(reecb.y().getData()[index(Var::PMEAS)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::QREF)], static_cast<ScalarT>(0.1), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::PORD)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(pe_node.read(), static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(qgen_node.read(), static_cast<ScalarT>(0.05), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.05), kTol);
        success *= isEqual(reecb.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.25), kTol);
        success *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome zeroTime()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeData();
        data.parameters[Params::Trv] = static_cast<RealT>(0.0);
        data.parameters[Params::Tp]  = static_cast<RealT>(0.0);

        ReecbT reecb(&bus, data);
        success                                *= (reecb.allocate() == 0);
        success                                *= (reecb.verify() == 0);
        reecb.y().getData()[index(Var::IQCMD)]  = static_cast<ScalarT>(0.2);
        reecb.y().getData()[index(Var::IPCMD)]  = static_cast<ScalarT>(0.6);
        reecb.y().setDataUpdated();
        success *= (reecb.initialize() == 0);
        success *= (reecb.tagDifferentiable() == 0);
        success *= (reecb.tag()[index(Var::VMEAS)] == true);
        success *= (reecb.tag()[index(Var::PMEAS)] == true);

        reecb.yp().getData()[index(Var::VMEAS)] = static_cast<ScalarT>(1.0);
        reecb.yp().getData()[index(Var::PMEAS)] = static_cast<ScalarT>(2.0);
        reecb.yp().setDataUpdated();
        success *= (reecb.evaluateResidual() == 0);
        success *= isEqual(reecb.getResidual().getData()[index(Var::VMEAS)], static_cast<ScalarT>(-1.0), kTol);
        success *= isEqual(reecb.getResidual().getData()[index(Var::PMEAS)], static_cast<ScalarT>(-2.0), kTol);

        reecb.yp().getData()[index(Var::VMEAS)] = ZERO<RealT>;
        reecb.y().getData()[index(Var::VMEAS)]  = static_cast<ScalarT>(0.99);
        reecb.y().setDataUpdated();
        reecb.yp().setDataUpdated();
        success *= (reecb.evaluateResidual() == 0);
        success *= isEqual(reecb.getResidual().getData()[index(Var::VMEAS)], static_cast<ScalarT>(10.0), kTol);

        return success.report(__func__);
      }

      TestOutcome qPriority()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                       = makeData();
        data.parameters[Params::Pqflag] = static_cast<IdxT>(0);
        data.parameters[Params::Imax]   = static_cast<RealT>(1.1);

        ReecbT reecb(&bus, data);
        success                                *= (reecb.allocate() == 0);
        success                                *= (reecb.verify() == 0);
        reecb.y().getData()[index(Var::IQCMD)]  = static_cast<ScalarT>(0.3);
        reecb.y().getData()[index(Var::IPCMD)]  = static_cast<ScalarT>(0.9);
        reecb.y().setDataUpdated();
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        const ScalarT ipmax  = std::sqrt(static_cast<ScalarT>(1.1 * 1.1 - 0.3 * 0.3));
        success             *= isEqual(reecb.y().getData()[index(Var::IQMAX)], static_cast<ScalarT>(1.1), kTol);
        success             *= isEqual(reecb.y().getData()[index(Var::IPMAX)], ipmax, kTol);
        success             *= isEqual(reecb.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.3), kTol);
        success             *= isEqual(reecb.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.9), kTol);
        success             *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome pPriority()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                       = makeData();
        data.parameters[Params::Pqflag] = static_cast<IdxT>(1);
        data.parameters[Params::Imax]   = static_cast<RealT>(1.1);

        ReecbT reecb(&bus, data);
        success                                *= (reecb.allocate() == 0);
        success                                *= (reecb.verify() == 0);
        reecb.y().getData()[index(Var::IQCMD)]  = static_cast<ScalarT>(0.3);
        reecb.y().getData()[index(Var::IPCMD)]  = static_cast<ScalarT>(0.9);
        reecb.y().setDataUpdated();
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        const ScalarT iqmax  = std::sqrt(static_cast<ScalarT>(1.1 * 1.1 - 0.9 * 0.9));
        success             *= isEqual(reecb.y().getData()[index(Var::IPMAX)], static_cast<ScalarT>(1.1), kTol);
        success             *= isEqual(reecb.y().getData()[index(Var::IQMAX)], iqmax, kTol);
        success             *= isEqual(reecb.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.3), kTol);
        success             *= isEqual(reecb.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.9), kTol);
        success             *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome voltageBand()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                      = makeData();
        data.parameters[Params::QFlag] = static_cast<IdxT>(1);
        data.parameters[Params::Vref0] = static_cast<RealT>(1.0);
        data.parameters[Params::Vdip]  = static_cast<RealT>(0.9);
        data.parameters[Params::Vup]   = static_cast<RealT>(1.1);
        data.parameters[Params::dbd1]  = static_cast<RealT>(0.0);
        data.parameters[Params::dbd2]  = static_cast<RealT>(0.0);
        data.parameters[Params::kqv]   = static_cast<RealT>(1.0);
        data.parameters[Params::Imax]  = static_cast<RealT>(2.0);

        ReecbT reecb(&bus, data);
        success                                *= (reecb.allocate() == 0);
        success                                *= (reecb.verify() == 0);
        reecb.y().getData()[index(Var::IQCMD)]  = ZERO<RealT>;
        reecb.y().getData()[index(Var::IPCMD)]  = static_cast<ScalarT>(0.5);
        reecb.y().setDataUpdated();
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        const ScalarT expected_sdip   = Math::inside(reecb.y().getData()[index(Var::VT)], static_cast<RealT>(0.9), static_cast<RealT>(1.1));
        const ScalarT expected_iqraw  = reecb.y().getData()[index(Var::IQBASE)] + (ONE<RealT> - reecb.y().getData()[index(Var::SDIP)]) * reecb.y().getData()[index(Var::IQV)];
        success                      *= isEqual(reecb.y().getData()[index(Var::SDIP)], expected_sdip, kTol);
        success                      *= isEqual(reecb.y().getData()[index(Var::IQRAW)], expected_iqraw, kTol);
        success                      *= (reecb.y().getData()[index(Var::SDIP)] < static_cast<ScalarT>(0.01));
        success                      *= (reecb.y().getData()[index(Var::IQV)] > static_cast<ScalarT>(0.1));
        success                      *= allZero(reecb);

        return success.report(__func__);
      }

      TestOutcome piSaturation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.13, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                       = makeData();
        data.parameters[Params::QFlag]  = static_cast<IdxT>(0);
        data.parameters[Params::Pqflag] = static_cast<IdxT>(0);
        data.parameters[Params::Vmin]   = static_cast<RealT>(0.9);
        data.parameters[Params::Vmax]   = static_cast<RealT>(1.05);
        data.parameters[Params::Kvp]    = static_cast<RealT>(10.0);
        data.parameters[Params::Kvi]    = static_cast<RealT>(60.0);
        data.parameters[Params::Vup]    = static_cast<RealT>(99.0);
        data.parameters[Params::Vdip]   = static_cast<RealT>(-99.0);
        data.parameters[Params::Imax]   = static_cast<RealT>(1.1);

        ReecbT reecb(&bus, data);
        success                                *= (reecb.allocate() == 0);
        success                                *= (reecb.verify() == 0);
        reecb.y().getData()[index(Var::IQCMD)]  = static_cast<ScalarT>(0.15 / 1.13);
        reecb.y().getData()[index(Var::IPCMD)]  = static_cast<ScalarT>(0.5 / 1.13);
        reecb.y().setDataUpdated();
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        const ScalarT piv_arg  = static_cast<ScalarT>(10.0) * reecb.y().getData()[index(Var::EPIV)] + reecb.y().getData()[index(Var::XPIV)];
        success               *= (piv_arg < -reecb.y().getData()[index(Var::IQMAX)]);
        success               *= allZero(reecb);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ReecbT reecb(&bus, makeData());
        success                                *= (reecb.allocate() == 0);
        success                                *= (reecb.verify() == 0);
        reecb.y().getData()[index(Var::IQCMD)]  = static_cast<ScalarT>(0.2);
        reecb.y().getData()[index(Var::IPCMD)]  = static_cast<ScalarT>(0.6);
        reecb.y().setDataUpdated();
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);
        success *= (reecb.evaluateJacobian() == 0);

        auto* jac  = reecb.getCooJacobian();
        success   *= (jac != nullptr);
        if (jac != nullptr)
        {
          success            *= (jac->getNnz() > 0);
          const auto* values  = jac->getValues();
          for (IdxT i = 0; i < jac->getNnz(); ++i)
          {
            success *= std::isfinite(values[i]);
          }
        }

        return success.report(__func__);
      }
#endif

      TestOutcome json()
      {
        TestStatus success = true;

        std::istringstream input(R"json(
{
  "header": {
    "format_version": 0,
    "format_revision": 1,
    "case_name": "renewable electrical control",
    "case_description": "REECB parser test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "params": { "kv": 1.0 } }
  ],
  "signals": [
    { "signal_id": 12, "name": "Iqcmd" },
    { "signal_id": 13, "name": "Ipcmd" }
  ],
  "devices": [
    {
      "class": "Reecb",
      "ports": { "bus": 1, "iqcmd": 12, "ipcmd": 13 },
      "id": "REE1",
      "params": {
        "mva": 100.0, "PfFlag": 0, "VFlag": 1, "QFlag": 0, "Pqflag": 1,
        "Trv": 0.0, "Tp": 0.02, "Vdip": 0.7, "Vup": 1.2,
        "dbd1": -0.01, "dbd2": 0.01, "kqv": 0.0, "Iql1": -1.0, "Iqh1": 1.0,
        "Qmax": 1.0, "Qmin": -1.0, "Kqp": 1.0, "Kqi": 0.0,
        "Vmax": 1.2, "Vmin": 0.8, "Kvp": 1.0, "Kvi": 0.0,
        "Tiq": 0.02, "Tpord": 0.02, "dPmax": 1.0, "dPmin": -1.0,
        "Pmax": 1.0, "Pmin": 0.0, "Imax": 2.0
      }
    }
  ]
}
)json");

        auto data  = PhasorDynamics::parseSystemModelData(input);
        success   *= (data.reecb.size() == 1);
        success   *= (std::get<IdxT>(data.reecb[0].parameters.at(Params::Pqflag)) == static_cast<IdxT>(1));
        success   *= (std::get<RealT>(data.reecb[0].parameters.at(Params::mva)) == static_cast<RealT>(100.0));
        success   *= (data.reecb[0].buses.at(Buses::bus) == static_cast<IdxT>(1));
        success   *= data.reecb[0].signal_inputs.empty();
        success   *= (data.reecb[0].signal_outputs.at(Outputs::iqcmd) == static_cast<IdxT>(12));
        success   *= (data.reecb[0].signal_outputs.at(Outputs::ipcmd) == static_cast<IdxT>(13));

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() == 27);

        return success.report(__func__);
      }

    private:
      static size_t index(Var variable)
      {
        return static_cast<size_t>(variable);
      }

      TestStatus allZero(const ReecbT& reecb) const
      {
        TestStatus success = true;

        for (size_t i = 0; i < reecb.getResidual().getSize(); ++i)
        {
          success *= isEqual(reecb.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reecb.yp().getData()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success;
      }

      auto makeMinimalData() -> PhasorDynamics::Converter::ReecbData<RealT, IdxT>
      {
        using Mon = PhasorDynamics::Converter::ReecbMonitorableVariables;

        PhasorDynamics::Converter::ReecbData<RealT, IdxT> data;
        data.device_class          = "Reecb";
        data.disambiguation_string = "reecb_test";
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);
        return data;
      }

      auto makeData() -> PhasorDynamics::Converter::ReecbData<RealT, IdxT>
      {
        auto data = makeMinimalData();

        data.parameters[Params::mva]    = static_cast<RealT>(100.0);
        data.parameters[Params::PfFlag] = static_cast<IdxT>(0);
        data.parameters[Params::VFlag]  = static_cast<IdxT>(1);
        data.parameters[Params::QFlag]  = static_cast<IdxT>(0);
        data.parameters[Params::Pqflag] = static_cast<IdxT>(1);
        data.parameters[Params::Trv]    = static_cast<RealT>(0.0);
        data.parameters[Params::Tp]     = static_cast<RealT>(0.02);
        data.parameters[Params::Vdip]   = static_cast<RealT>(0.7);
        data.parameters[Params::Vup]    = static_cast<RealT>(1.2);
        data.parameters[Params::dbd1]   = static_cast<RealT>(-0.01);
        data.parameters[Params::dbd2]   = static_cast<RealT>(0.01);
        data.parameters[Params::kqv]    = static_cast<RealT>(0.0);
        data.parameters[Params::Iql1]   = static_cast<RealT>(-1.0);
        data.parameters[Params::Iqh1]   = static_cast<RealT>(1.0);
        data.parameters[Params::Qmax]   = static_cast<RealT>(1.0);
        data.parameters[Params::Qmin]   = static_cast<RealT>(-1.0);
        data.parameters[Params::Kqp]    = static_cast<RealT>(1.0);
        data.parameters[Params::Kqi]    = static_cast<RealT>(0.0);
        data.parameters[Params::Vmax]   = static_cast<RealT>(1.2);
        data.parameters[Params::Vmin]   = static_cast<RealT>(0.8);
        data.parameters[Params::Kvp]    = static_cast<RealT>(1.0);
        data.parameters[Params::Kvi]    = static_cast<RealT>(0.0);
        data.parameters[Params::Tiq]    = static_cast<RealT>(0.02);
        data.parameters[Params::Tpord]  = static_cast<RealT>(0.02);
        data.parameters[Params::dPmax]  = static_cast<RealT>(1.0);
        data.parameters[Params::dPmin]  = static_cast<RealT>(-1.0);
        data.parameters[Params::Pmax]   = static_cast<RealT>(1.0);
        data.parameters[Params::Pmin]   = static_cast<RealT>(0.0);
        data.parameters[Params::Imax]   = static_cast<RealT>(2.0);

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
