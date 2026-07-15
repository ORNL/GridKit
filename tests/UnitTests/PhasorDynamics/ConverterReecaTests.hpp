#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/Reeca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECA/ReecaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ConverterReecaTests
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

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeReecaData());
        success *=
            (reeca.size()
             == static_cast<IdxT>(PhasorDynamics::Converter::ReecaInternalVariables::MAXIMUM));
        success *= (reeca.getMonitor() != nullptr);
        success *= (reeca.verify() == 0);

        auto negative_vdip_reeca = makeReecaData();
        negative_vdip_reeca.parameters[PhasorDynamics::Converter::ReecaParameters::Vdip] =
            static_cast<RealT>(-0.1);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> negative_vdip_reeca_model(
            &bus,
            negative_vdip_reeca);
        success *= (negative_vdip_reeca_model.verify() > 0);

        auto bad_reeca_band = makeReecaData();
        bad_reeca_band.parameters[PhasorDynamics::Converter::ReecaParameters::Vdip] =
            static_cast<RealT>(1.2);
        bad_reeca_band.parameters[PhasorDynamics::Converter::ReecaParameters::Vup] =
            static_cast<RealT>(1.2);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_reeca_band_model(
            &bus,
            bad_reeca_band);
        success *= (bad_reeca_band_model.verify() > 0);

        auto bad_reeca = makeReecaData();
        bad_reeca.parameters[PhasorDynamics::Converter::ReecaParameters::Imax] =
            static_cast<RealT>(-1.0);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_reeca_model(&bus, bad_reeca);
        success *= (bad_reeca_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> defaults_model(
            &bus,
            makeMinimalReecaData());
        success *= (defaults_model.verify() == 0);
        success *= (bus.allocate() == 0);
        success *= (bus.initialize() == 0);
        success *= (defaults_model.allocate() == 0);
        success *= (defaults_model.initialize() == 0);
        success *= (defaults_model.evaluateResidual() == 0);
        for (IdxT i = 0; i < defaults_model.getResidual().getSize(); ++i)
        {
          success *= isEqual(
              defaults_model.getResidual().getData()[static_cast<size_t>(i)],
              static_cast<ScalarT>(0.0),
              kTol);
        }

        auto missing_mva = makeMinimalReecaData();
        missing_mva.parameters.erase(Params::mva);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> missing_mva_model(
            &bus,
            missing_mva);
        success *= (missing_mva_model.verify() > 0);

        auto bad_switch                       = makeReecaData();
        bad_switch.parameters[Params::PfFlag] = static_cast<IdxT>(2);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_switch_model(&bus, bad_switch);
        success *= (bad_switch_model.verify() > 0);

        auto bad_pflag                      = makeReecaData();
        bad_pflag.parameters[Params::PFlag] = static_cast<IdxT>(2);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> bad_pflag_model(&bus, bad_pflag);
        success *= (bad_pflag_model.verify() > 0);

        auto missing_speed                      = makeReecaData();
        missing_speed.parameters[Params::PFlag] = static_cast<IdxT>(1);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> missing_speed_model(
            &bus,
            missing_speed);
        success *= (missing_speed_model.verify() > 0);

        success *= invalidParameterCase(bus, Params::mva, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Trv, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Tp, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Vdip, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Vup, static_cast<RealT>(0.5));
        success *= invalidParameterCase(bus, Params::dbd1, static_cast<RealT>(0.1));
        success *= invalidParameterCase(bus, Params::dbd2, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Iql1, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Thld, static_cast<RealT>(1.0));
        success *= invalidParameterCase(bus, Params::Thld2, static_cast<RealT>(1.0));
        success *= invalidParameterCase(bus, Params::Qmin, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Vmin, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Tiq, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Tpord, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::dPmin, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::dPmax, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Pmin, static_cast<RealT>(2.0));
        success *= invalidParameterCase(bus, Params::Imax, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Vq2, static_cast<RealT>(0.2));
        success *= invalidParameterCase(bus, Params::Iq1, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::Vp2, static_cast<RealT>(0.2));
        success *= invalidParameterCase(bus, Params::Ip1, static_cast<RealT>(-0.1));

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> valid_model(&bus, makeReecaData());
        success *= (valid_model.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome reecaSignalsInitializationAndResidual()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;
        using Ext = PhasorDynamics::Converter::ReecaExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        ScalarT pe_value{0.75};
        ScalarT qgen_value{0.2};
        ScalarT omega_value{1.01};
        ScalarT pref_value{99.0};
        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.75};
        IdxT    pe_index    = 20;
        IdxT    qgen_index  = 21;
        IdxT    omega_index = 22;
        IdxT    pref_index  = 23;
        IdxT    iqcmd_index = 24;
        IdxT    ipcmd_index = 25;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        omega_node.set(&omega_value, &omega_index);
        pref_node.set(&pref_value, &pref_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        auto data                                                          = makeReecaData();
        data.parameters[PhasorDynamics::Converter::ReecaParameters::QFlag] = static_cast<IdxT>(1);
        data.parameters[PhasorDynamics::Converter::ReecaParameters::PFlag] = static_cast<IdxT>(1);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reeca.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        reeca.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(qgen_value);
        ipcmd_node.init(pe_value);
        success *= (reeca.initialize() == 0);
        success *= (reeca.tagDifferentiable() == 0);
        success *= (reeca.evaluateResidual() == 0);

        success                     *= isEqual(reeca.y().getData()[index(Var::VMEAS)], static_cast<ScalarT>(1.0), kTol);
        success                     *= isEqual(reeca.y().getData()[index(Var::PMEAS)], pe_value, kTol);
        success                     *= isEqual(reeca.y().getData()[index(Var::QREF)], qgen_value, kTol);
        success                     *= isEqual(reeca.y().getData()[index(Var::PORD)], pe_value, kTol);
        const ScalarT expected_pref  = pe_value / static_cast<ScalarT>(1.01);
        success                     *= isEqual(pref_node.read(), expected_pref, kTol);
        success                     *= isEqual(iqcmd_node.read(), reeca.y().getData()[index(Var::IQCMD)], kTol);
        success                     *= isEqual(ipcmd_node.read(), reeca.y().getData()[index(Var::IPCMD)], kTol);
        success                     *= (reeca.tag()[index(Var::VMEAS)] == true);
        success                     *= (reeca.tag()[index(Var::PMEAS)] == true);

        for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
        {
          success *= isEqual(reeca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reeca.yp().getData()[i], static_cast<ScalarT>(0.0), kTol);
        }

        omega_value  = static_cast<ScalarT>(1.02);
        success     *= (reeca.evaluateResidual() == 0);
        const ScalarT expected_fpord =
            (omega_value * expected_pref - reeca.y().getData()[index(Var::PORD)])
            / static_cast<RealT>(0.02);
        success *= isEqual(
            reeca.getResidual().getData()[index(Var::FPORD)],
            expected_fpord,
            kTol);

        return success.report(__func__);
      }

      TestOutcome rejectsHalfConnectedElectricalFeedback()
      {
        using Ext = PhasorDynamics::Converter::ReecaExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        ScalarT signal_value{0.6};
        IdxT    signal_index = 24;

        PhasorDynamics::SignalNode<ScalarT, IdxT> signal_node;
        signal_node.set(&signal_value, &signal_index);

        {
          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeReecaData());
          reeca.getSignals().template attachSignalNode<Ext::PE>(&signal_node);
          success *= (reeca.verify() > 0);
        }

        {
          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeReecaData());
          reeca.getSignals().template attachSignalNode<Ext::QGEN>(&signal_node);
          success *= (reeca.verify() > 0);
        }

        {
          PhasorDynamics::SignalNode<ScalarT, IdxT>       unlinked_node;
          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeReecaData());
          reeca.getSignals().template attachSignalNode<Ext::OMEGA>(&unlinked_node);
          success *= (reeca.verify() > 0);
        }

        return success.report(__func__);
      }

      TestOutcome reecaCommandSignalInitialization()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.6};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeReecaData());
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.6));
        success *= (reeca.initialize() == 0);
        success *= (reeca.tagDifferentiable() == 0);
        success *= (reeca.evaluateResidual() == 0);

        success *= isEqual(reeca.y().getData()[index(Var::PMEAS)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::QREF)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::PORD)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.2), kTol);

        for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
        {
          success *= isEqual(reeca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reeca.yp().getData()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome reecaElectricalFeedbackUsesMvaBase()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;
        using Ext = PhasorDynamics::Converter::ReecaExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeReecaData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        ScalarT pe_value{0.25};
        ScalarT qgen_value{0.05};
        ScalarT qext_value{99.0};
        ScalarT pfaref_value{99.0};
        ScalarT pref_value{99.0};
        ScalarT iqcmd_value{0.0};
        ScalarT ipcmd_value{0.0};
        IdxT    pe_index     = 25;
        IdxT    qgen_index   = 26;
        IdxT    qext_index   = 27;
        IdxT    pfaref_index = 28;
        IdxT    pref_index   = 29;
        IdxT    iqcmd_index  = 30;
        IdxT    ipcmd_index  = 31;

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

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reeca.getSignals().template attachSignalNode<Ext::QEXT>(&qext_node);
        reeca.getSignals().template attachSignalNode<Ext::PFAREF>(&pfaref_node);
        reeca.getSignals().template attachSignalNode<Ext::PREF>(&pref_node);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.05));
        ipcmd_node.init(static_cast<ScalarT>(0.25));
        success *= (reeca.initialize() == 0);
        success *= (reeca.evaluateResidual() == 0);

        const ScalarT expected_pfaref =
            static_cast<ScalarT>(std::atan(static_cast<RealT>(0.1 / 0.5)));
        success *= isEqual(reeca.y().getData()[index(Var::PMEAS)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::QREF)], static_cast<ScalarT>(0.1), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::PORD)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(pe_node.read(), static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(qgen_node.read(), static_cast<ScalarT>(0.05), kTol);
        success *= isEqual(qext_node.read(), static_cast<ScalarT>(0.05), kTol);
        success *= isEqual(pfaref_node.read(), expected_pfaref, kTol);
        success *= isEqual(pref_node.read(), static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(iqcmd_node.read(), static_cast<ScalarT>(0.05), kTol);
        success *= isEqual(ipcmd_node.read(), static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.25), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.05), kTol);

        for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
        {
          success *= isEqual(reeca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome reecaReferenceFallbackAtAngle()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        auto data                       = makeReecaData();
        data.parameters[Params::PfFlag] = static_cast<IdxT>(1);
        data.parameters[Params::Qmin]   = static_cast<RealT>(-2.0);
        data.parameters[Params::Qmax]   = static_cast<RealT>(2.0);
        data.parameters[Params::Pmax]   = static_cast<RealT>(2.0);
        data.parameters[Params::Imax]   = static_cast<RealT>(2.0);

        ScalarT iqcmd_value{-0.2};
        ScalarT ipcmd_value{1.6};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(-0.2));
        ipcmd_node.init(static_cast<ScalarT>(1.6));
        success *= (reeca.initialize() == 0);
        success *= (reeca.tagDifferentiable() == 0);
        success *= (reeca.evaluateResidual() == 0);

        success *= isEqual(reeca.y().getData()[index(Var::PORD)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::QREF)], static_cast<ScalarT>(-0.2), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(-0.2), kTol);

        for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
        {
          success *= isEqual(reeca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reeca.yp().getData()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome zeroTimeConstants()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeReecaData();
        data.parameters[Params::Trv] = static_cast<RealT>(0.0);
        data.parameters[Params::Tp]  = static_cast<RealT>(0.0);

        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.6};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        success *= (reeca.tagDifferentiable() == 0);
        success *= (reeca.tag()[index(Var::VMEAS)]);
        success *= (reeca.tag()[index(Var::PMEAS)]);
        success *= (!reeca.tag()[index(Var::VT)]);
        success *= (reeca.tag()[index(Var::XPIQ)]);
        success *= (reeca.tag()[index(Var::XPIV)]);
        success *= (reeca.tag()[index(Var::QV)]);
        success *= (reeca.tag()[index(Var::PORD)]);

        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.6));
        success *= (reeca.initialize() == 0);

        reeca.yp().getData()[index(Var::VMEAS)] = static_cast<ScalarT>(1.0);
        reeca.yp().getData()[index(Var::PMEAS)] = static_cast<ScalarT>(2.0);
        reeca.yp().setDataUpdated();
        success *= (reeca.evaluateResidual() == 0);
        success *= isEqual(
            reeca.getResidual().getData()[index(Var::VMEAS)],
            static_cast<ScalarT>(-1.0),
            kTol);
        success *= isEqual(
            reeca.getResidual().getData()[index(Var::PMEAS)],
            static_cast<ScalarT>(-2.0),
            kTol);

        reeca.y().getData()[index(Var::VMEAS)]  = static_cast<ScalarT>(0.99);
        reeca.y().getData()[index(Var::PMEAS)]  = static_cast<ScalarT>(0.59);
        reeca.yp().getData()[index(Var::VMEAS)] = ZERO<RealT>;
        reeca.yp().getData()[index(Var::PMEAS)] = ZERO<RealT>;
        reeca.y().setDataUpdated();
        reeca.yp().setDataUpdated();
        success *= (reeca.evaluateResidual() == 0);
        success *= isEqual(
            reeca.getResidual().getData()[index(Var::VMEAS)],
            static_cast<ScalarT>(10.0),
            kTol);
        success *= isEqual(
            reeca.getResidual().getData()[index(Var::PMEAS)],
            static_cast<ScalarT>(10.0),
            kTol);

        return success.report(__func__);
      }

      TestOutcome nonnegativeCurrentLimits()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto runCase = [&](const IdxT  priority,
                           const RealT imax,
                           const RealT iq_limit,
                           const RealT ip_limit)
        {
          auto data                       = makeReecaData();
          data.parameters[Params::Pqflag] = priority;
          data.parameters[Params::Imax]   = imax;
          data.parameters[Params::Iq1]    = iq_limit;
          data.parameters[Params::Iq2]    = iq_limit;
          data.parameters[Params::Iq3]    = iq_limit;
          data.parameters[Params::Iq4]    = iq_limit;
          data.parameters[Params::Ip1]    = ip_limit;
          data.parameters[Params::Ip2]    = ip_limit;
          data.parameters[Params::Ip3]    = ip_limit;
          data.parameters[Params::Ip4]    = ip_limit;

          ScalarT iqcmd_value{ZERO<RealT>};
          ScalarT ipcmd_value{ZERO<RealT>};
          IdxT    iqcmd_index = 22;
          IdxT    ipcmd_index = 23;

          PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
          iqcmd_node.set(&iqcmd_value, &iqcmd_index);
          ipcmd_node.set(&ipcmd_value, &ipcmd_index);

          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
          reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
          reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

          success *= (reeca.allocate() == 0);
          success *= (reeca.verify() == 0);
          iqcmd_node.init(ZERO<RealT>);
          ipcmd_node.init(ZERO<RealT>);
          success *= (reeca.initialize() == 0);
          success *= (reeca.evaluateResidual() == 0);

          const auto*   y = reeca.y().getData();
          const ScalarT expected_iq_limit =
              static_cast<ScalarT>(iq_limit) * y[index(Var::IQCIRC)]
              / Math::max(static_cast<ScalarT>(iq_limit), y[index(Var::IQCIRC)]);
          const ScalarT expected_ip_limit =
              static_cast<ScalarT>(ip_limit) * y[index(Var::IPCIRC)]
              / Math::max(static_cast<ScalarT>(ip_limit), y[index(Var::IPCIRC)]);
          const ScalarT iq_circle = y[index(Var::IQCIRC)] * y[index(Var::IQCIRC)]
                                    + static_cast<ScalarT>(priority)
                                          * y[index(Var::IPCMD)] * y[index(Var::IPCMD)];
          const ScalarT ip_circle = y[index(Var::IPCIRC)] * y[index(Var::IPCIRC)]
                                    + (ONE<RealT> - static_cast<ScalarT>(priority))
                                          * y[index(Var::IQCMD)] * y[index(Var::IQCMD)];
          success *= isEqual(iq_circle, static_cast<ScalarT>(imax * imax), kTol);
          success *= isEqual(ip_circle, static_cast<ScalarT>(imax * imax), kTol);
          success *= isEqual(y[index(Var::IQMAX)], expected_iq_limit, kTol);
          success *= isEqual(y[index(Var::IPMAX)], expected_ip_limit, kTol);
          success *= (y[index(Var::IQMAX)] >= ZERO<RealT>);
          success *= (y[index(Var::IPMAX)] >= ZERO<RealT>);
          success *= (y[index(Var::IQCMD)] >= -y[index(Var::IQMAX)]);
          success *= (y[index(Var::IQCMD)] <= y[index(Var::IQMAX)]);
          success *= (y[index(Var::IPCMD)] >= ZERO<RealT>);
          success *= (y[index(Var::IPCMD)] <= y[index(Var::IPMAX)]);
          for (IdxT i = 0; i < reeca.getResidual().getSize(); ++i)
          {
            success *= isEqual(
                reeca.getResidual().getData()[static_cast<size_t>(i)],
                ZERO<RealT>,
                kTol);
          }
        };

        for (IdxT priority : {static_cast<IdxT>(0), static_cast<IdxT>(1)})
        {
          runCase(priority, ZERO<RealT>, ZERO<RealT>, ZERO<RealT>);
          runCase(
              priority,
              ZERO<RealT>,
              static_cast<RealT>(2.0),
              static_cast<RealT>(2.0));
          runCase(
              priority,
              static_cast<RealT>(2.0),
              ZERO<RealT>,
              static_cast<RealT>(2.0));
          runCase(
              priority,
              static_cast<RealT>(2.0),
              static_cast<RealT>(2.0),
              ZERO<RealT>);
          runCase(
              priority,
              static_cast<RealT>(1.0e-3),
              static_cast<RealT>(1.0e-3),
              static_cast<RealT>(1.0e-3));
        }

        return success.report(__func__);
      }

      TestOutcome voltageDipAndVdl()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                     = makeReecaData();
        data.parameters[Params::Vdip] = static_cast<RealT>(0.9);
        data.parameters[Params::Vup]  = static_cast<RealT>(1.1);
        data.parameters[Params::kqv]  = static_cast<RealT>(2.0);
        data.parameters[Params::Iql1] = static_cast<RealT>(-2.0);
        data.parameters[Params::Iqh1] = static_cast<RealT>(2.0);
        data.parameters[Params::Iq1]  = static_cast<RealT>(1.2);
        data.parameters[Params::Iq2]  = static_cast<RealT>(1.1);
        data.parameters[Params::Iq3]  = static_cast<RealT>(1.0);
        data.parameters[Params::Iq4]  = static_cast<RealT>(0.9);
        data.parameters[Params::Ip1]  = static_cast<RealT>(1.3);
        data.parameters[Params::Ip2]  = static_cast<RealT>(1.2);
        data.parameters[Params::Ip3]  = static_cast<RealT>(1.1);
        data.parameters[Params::Ip4]  = static_cast<RealT>(1.0);

        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.6};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.6));
        success *= (reeca.initialize() == 0);

        bus.Vr() = static_cast<ScalarT>(0.8);
        bus.Vi() = ZERO<RealT>;
        bus.y().setDataUpdated();

        auto* y                     = reeca.y().getData();
        y[index(Var::VMEAS)]        = static_cast<ScalarT>(0.8);
        y[index(Var::VT)]           = static_cast<ScalarT>(0.8);
        y[index(Var::VMEASSAFE)]    = static_cast<ScalarT>(0.8);
        const ScalarT expected_sdip = Math::inside(
            y[index(Var::VT)],
            static_cast<RealT>(0.9),
            static_cast<RealT>(1.1));
        const ScalarT expected_verr = Math::deadband2(
            static_cast<ScalarT>(1.0) - y[index(Var::VMEAS)],
            static_cast<RealT>(-0.01),
            static_cast<RealT>(0.01));
        const ScalarT expected_iqv = Math::clamp(
            static_cast<ScalarT>(2.0) * expected_verr,
            static_cast<RealT>(-2.0),
            static_cast<RealT>(2.0));
        y[index(Var::SDIP)] = expected_sdip;
        y[index(Var::VERR)] = expected_verr;
        y[index(Var::IQV)]  = expected_iqv;
        reeca.y().setDataUpdated();

        const ScalarT expected_gq = static_cast<ScalarT>(1.2)
                                    + Math::linseg(
                                        y[index(Var::VMEAS)],
                                        static_cast<RealT>(0.2),
                                        static_cast<RealT>(0.5),
                                        static_cast<RealT>(-0.1))
                                    + Math::linseg(
                                        y[index(Var::VMEAS)],
                                        static_cast<RealT>(0.5),
                                        static_cast<RealT>(0.75),
                                        static_cast<RealT>(-0.1))
                                    + Math::linseg(
                                        y[index(Var::VMEAS)],
                                        static_cast<RealT>(0.75),
                                        static_cast<RealT>(1.0),
                                        static_cast<RealT>(-0.1));
        const ScalarT expected_gp = static_cast<ScalarT>(1.3)
                                    + Math::linseg(
                                        y[index(Var::VMEAS)],
                                        static_cast<RealT>(0.2),
                                        static_cast<RealT>(0.5),
                                        static_cast<RealT>(-0.1))
                                    + Math::linseg(
                                        y[index(Var::VMEAS)],
                                        static_cast<RealT>(0.5),
                                        static_cast<RealT>(0.75),
                                        static_cast<RealT>(-0.1))
                                    + Math::linseg(
                                        y[index(Var::VMEAS)],
                                        static_cast<RealT>(0.75),
                                        static_cast<RealT>(1.0),
                                        static_cast<RealT>(-0.1));
        const ScalarT expected_iqmax_residual =
            -y[index(Var::IQMAX)]
            + expected_gq * y[index(Var::IQCIRC)]
                  / Math::max(expected_gq, y[index(Var::IQCIRC)]);
        const ScalarT expected_ipmax_residual =
            -y[index(Var::IPMAX)]
            + expected_gp * y[index(Var::IPCIRC)]
                  / Math::max(expected_gp, y[index(Var::IPCIRC)]);
        const ScalarT expected_iqraw_residual =
            -y[index(Var::IQRAW)] + y[index(Var::QV)]
            + (ONE<RealT> - expected_sdip) * expected_iqv;

        success *= (reeca.evaluateResidual() == 0);

        const auto* f  = reeca.getResidual().getData();
        success       *= (std::abs(static_cast<RealT>(expected_sdip)) < kTol);
        success       *= (static_cast<RealT>(expected_iqv) > static_cast<RealT>(0.1));
        success       *= (std::abs(static_cast<RealT>(expected_gq) - static_cast<RealT>(0.98))
                    < static_cast<RealT>(1.0e-6));
        success       *= (std::abs(static_cast<RealT>(expected_gp) - static_cast<RealT>(1.08))
                    < static_cast<RealT>(1.0e-6));
        success       *= isEqual(f[index(Var::SDIP)], ZERO<RealT>, kTol);
        success       *= isEqual(f[index(Var::VERR)], ZERO<RealT>, kTol);
        success       *= isEqual(f[index(Var::IQV)], ZERO<RealT>, kTol);
        success       *= isEqual(f[index(Var::IQMAX)], expected_iqmax_residual, kTol);
        success       *= isEqual(f[index(Var::IPMAX)], expected_ipmax_residual, kTol);
        success       *= isEqual(f[index(Var::IQRAW)], expected_iqraw_residual, kTol);
        success       *= isEqual(f[index(Var::XPIQ)], ZERO<RealT>, kTol);
        success       *= isEqual(f[index(Var::XPIV)], ZERO<RealT>, kTol);
        success       *= isEqual(f[index(Var::QV)], ZERO<RealT>, kTol);
        success       *= isEqual(f[index(Var::PORD)], ZERO<RealT>, kTol);

        return success.report(__func__);
      }

      TestOutcome outputAvailability()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        auto                               data = makeReecaData();
        addAllMonitors(data);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);

        success *= (reeca.getMonitor() != nullptr);
        success *= (!reeca.getMonitor()->empty());
        success *= (reeca.allocate() == 0);

        auto* y              = reeca.y().getData();
        y[index(Var::IQCMD)] = static_cast<ScalarT>(0.11);
        y[index(Var::IPCMD)] = static_cast<ScalarT>(0.22);
        y[index(Var::VMEAS)] = static_cast<ScalarT>(0.33);
        y[index(Var::PMEAS)] = static_cast<ScalarT>(0.44);
        reeca.y().setDataUpdated();

        RealT                                     time = ZERO<RealT>;
        Model::VariableMonitorController<ScalarT> controller(time);
        controller.addMonitor(reeca.getMonitor());

        std::stringstream os;
        controller.addSink({Model::VariableMonitorFormat::CSV}, os);
        controller.print();

        auto values  = Tokenizer<RealT>(os.str(), ',')();
        success     *= (values.size() == 5);
        if (values.size() == 5)
        {
          success *= isEqual(values[1], static_cast<RealT>(0.11), kTol);
          success *= isEqual(values[2], static_cast<RealT>(0.22), kTol);
          success *= isEqual(values[3], static_cast<RealT>(0.33), kTol);
          success *= isEqual(values[4], static_cast<RealT>(0.44), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome allocationAndAbsoluteTolerance()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>        bus(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, makeReecaData());
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);

        success                   *= (reeca.allocate() == 0);
        auto* y_data               = reeca.y().getData();
        y_data[index(Var::IQCMD)]  = static_cast<ScalarT>(0.37);
        reeca.y().setDataUpdated();

        const RealT absolute_tolerance  = static_cast<RealT>(2.5e-7);
        success                        *= (reeca.setAbsoluteTolerance(absolute_tolerance) == 0);
        success                        *= (reeca.absoluteTolerance().getSize() == reeca.size());
        for (IdxT i = 0; i < reeca.absoluteTolerance().getSize(); ++i)
        {
          success *= isEqual(
              reeca.absoluteTolerance().getData()[static_cast<size_t>(i)],
              static_cast<ScalarT>(absolute_tolerance),
              kTol);
        }

        success *= (reeca.allocate() == 0);
        success *= (reeca.y().getData() == y_data);
        success *= isEqual(
            reeca.y().getData()[index(Var::IQCMD)],
            static_cast<ScalarT>(0.37),
            kTol);
        success *= isEqual(iqcmd_node.read(), static_cast<ScalarT>(0.37), kTol);
        for (IdxT i = 0; i < reeca.absoluteTolerance().getSize(); ++i)
        {
          success *= isEqual(
              reeca.absoluteTolerance().getData()[static_cast<size_t>(i)],
              static_cast<ScalarT>(absolute_tolerance),
              kTol);
        }

        return success.report(__func__);
      }

      TestOutcome piInitializationContract()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        success *= (bus.allocate() == 0);
        success *= (bus.initialize() == 0);

        auto data                      = makeReecaData();
        data.parameters[Params::VFlag] = static_cast<IdxT>(0);
        data.parameters[Params::Vref1] = static_cast<RealT>(0.8000000005);
        data.parameters[Params::Vmin]  = static_cast<RealT>(0.0);
        data.parameters[Params::Vmax]  = static_cast<RealT>(2.0);
        data.parameters[Params::Kvi]   = static_cast<RealT>(1.0e10);
        data.parameters[Params::Imax]  = static_cast<RealT>(2.0);

        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.6};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(iqcmd_value);
        ipcmd_node.init(ipcmd_value);
        success *= (reeca.initialize() == 0);

        const auto*   y = reeca.y().getData();
        const ScalarT expected_vpiq =
            y[index(Var::QREF)] + static_cast<ScalarT>(0.8000000005);
        const ScalarT expected_epiv = expected_vpiq - y[index(Var::VMEAS)];
        const ScalarT xpiv_arg =
            y[index(Var::XPIV)] + y[index(Var::EPIV)];

        success *= isEqual(y[index(Var::VPIQ)], expected_vpiq, kTol);
        success *= isEqual(y[index(Var::EPIV)], expected_epiv, kTol);
        success *= isEqual(
            xpiv_arg,
            y[index(Var::IQMAX)] + static_cast<ScalarT>(0.1),
            kTol);

        success *= (reeca.evaluateResidual() == 0);
        for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
        {
          success *= isEqual(
              reeca.getResidual().getData()[i],
              static_cast<ScalarT>(0.0),
              kTol);
        }

        return success.report(__func__);
      }

      TestOutcome priorityInitialization()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                       = makeReecaData();
        data.parameters[Params::Pqflag] = static_cast<IdxT>(1);
        data.parameters[Params::Imax]   = static_cast<RealT>(1.2);

        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.75};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reeca.allocate() == 0);
        success *= (reeca.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.75));
        success *= (reeca.initialize() == 0);

        success *= isEqual(reeca.y().getData()[index(Var::VT)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IPCMD)], static_cast<ScalarT>(0.75), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IQCMD)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(reeca.y().getData()[index(Var::IPCIRC)], static_cast<ScalarT>(1.2), kTol);

        const auto circle =
            reeca.y().getData()[index(Var::IQCIRC)] * reeca.y().getData()[index(Var::IQCIRC)]
            + reeca.y().getData()[index(Var::IPCMD)] * reeca.y().getData()[index(Var::IPCMD)];
        success *= isEqual(circle, static_cast<ScalarT>(1.2 * 1.2), kTol);

        success *= (reeca.evaluateResidual() == 0);
        for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
        {
          success *= isEqual(reeca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome initializesSaturatedStartConsistently()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;

        TestStatus success = true;

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.13, 0.0);
          bus.allocate();
          bus.initialize();

          auto data                       = makeReecaData();
          data.parameters[Params::QFlag]  = static_cast<IdxT>(0);
          data.parameters[Params::VFlag]  = static_cast<IdxT>(0);
          data.parameters[Params::Pqflag] = static_cast<IdxT>(0);
          data.parameters[Params::Vmin]   = static_cast<RealT>(0.9);
          data.parameters[Params::Vmax]   = static_cast<RealT>(1.05);
          data.parameters[Params::Vref1]  = static_cast<RealT>(0.04);
          data.parameters[Params::Kvp]    = static_cast<RealT>(10.0);
          data.parameters[Params::Kvi]    = static_cast<RealT>(60.0);
          data.parameters[Params::Vup]    = static_cast<RealT>(99.0);
          data.parameters[Params::Vdip]   = static_cast<RealT>(0.0);
          data.parameters[Params::Imax]   = static_cast<RealT>(1.1);

          ScalarT iqcmd_value{0.15 / 1.13};
          ScalarT ipcmd_value{0.5 / 1.13};
          IdxT    iqcmd_index = 22;
          IdxT    ipcmd_index = 23;

          PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
          iqcmd_node.set(&iqcmd_value, &iqcmd_index);
          ipcmd_node.set(&ipcmd_value, &ipcmd_index);

          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
          reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
          reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

          success *= (reeca.allocate() == 0);
          success *= (reeca.verify() == 0);
          iqcmd_node.init(static_cast<ScalarT>(0.15 / 1.13));
          ipcmd_node.init(static_cast<ScalarT>(0.5 / 1.13));
          success *= (reeca.initialize() == 0);

          const ScalarT expected_epiv = reeca.y().getData()[index(Var::QREF)]
                                        + static_cast<ScalarT>(0.04)
                                        - reeca.y().getData()[index(Var::VMEAS)];
          success *= isEqual(reeca.y().getData()[index(Var::EPIV)], expected_epiv, kTol);

          const ScalarT piv_arg = static_cast<ScalarT>(10.0)
                                      * reeca.y().getData()[index(Var::EPIV)]
                                  + reeca.y().getData()[index(Var::XPIV)];
          success *= (piv_arg < -reeca.y().getData()[index(Var::IQMAX)]);

          success *= (reeca.evaluateResidual() == 0);
          for (size_t i = 0; i < reeca.getResidual().getSize(); ++i)
          {
            success *= isEqual(reeca.getResidual().getData()[i], static_cast<ScalarT>(0.0), kTol);
          }
        }

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.0);
          bus.allocate();
          bus.initialize();

          auto data                     = makeReecaData();
          data.parameters[Params::Vdip] = static_cast<RealT>(0.85);

          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
          success *= (reeca.allocate() == 0);
          success *= (reeca.verify() == 0);
          success *= (reeca.initialize() == 1);
        }

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
          bus.allocate();
          bus.initialize();

          auto data                     = makeReecaData();
          data.parameters[Params::Pmax] = static_cast<RealT>(2.0);
          data.parameters[Params::Ip1]  = static_cast<RealT>(1.0);
          data.parameters[Params::Ip2]  = static_cast<RealT>(1.0);
          data.parameters[Params::Ip3]  = static_cast<RealT>(1.0);
          data.parameters[Params::Ip4]  = static_cast<RealT>(1.0);

          ScalarT iqcmd_value{0.0};
          ScalarT ipcmd_value{1.5};
          IdxT    iqcmd_index = 22;
          IdxT    ipcmd_index = 23;

          PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
          iqcmd_node.set(&iqcmd_value, &iqcmd_index);
          ipcmd_node.set(&ipcmd_value, &ipcmd_index);

          PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);
          reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
          reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

          success *= (reeca.allocate() == 0);
          success *= (reeca.verify() == 0);
          iqcmd_node.init(static_cast<ScalarT>(0.0));
          ipcmd_node.init(static_cast<ScalarT>(1.5));
          success *= (reeca.initialize() == 1);
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto dependency_tracking_jacobian = DependencyTrackingJacobian();
        auto enzyme_jacobian              = EnzymeJacobian();

        success          *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
        const auto nrows  = std::min(dependency_tracking_jacobian.size(), enzyme_jacobian.size());
        for (size_t i = 0; i < nrows; ++i)
        {
          const auto dependency_tracking_row =
              pruneJacobianTails(dependency_tracking_jacobian[i], static_cast<RealT>(1.0e-8));
          const auto enzyme_row =
              pruneJacobianTails(enzyme_jacobian[i], static_cast<RealT>(1.0e-8));
          success *= isEqual(dependency_tracking_row, enzyme_row, static_cast<RealT>(1.0e-8));
        }

        return success.report(__func__);
      }
#endif

      TestOutcome jsonParseAndSystemAssembly()
      {
        TestStatus success = true;

        {
          std::istringstream input(commandOnlySystemJson());

          auto data  = PhasorDynamics::parseSystemModelData(input);
          success   *= (data.reeca.size() == 1);
          success   *= (std::get<IdxT>(
                          data.reeca[0].parameters.at(PhasorDynamics::Converter::ReecaParameters::Pqflag))
                      == static_cast<IdxT>(1));
          success   *= (std::get<RealT>(
                          data.reeca[0].parameters.at(PhasorDynamics::Converter::ReecaParameters::mva))
                      == static_cast<RealT>(100.0));
          success   *= (data.reeca[0].buses.at(Buses::bus) == static_cast<IdxT>(1));
          success   *= data.reeca[0].signal_inputs.empty();
          success   *= (data.reeca[0].signal_outputs.at(Outputs::iqcmd) == static_cast<IdxT>(12));
          success   *= (data.reeca[0].signal_outputs.at(Outputs::ipcmd) == static_cast<IdxT>(13));

          PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
          success *= (system.allocate() == 0);
          success *= (system.initialize() == 0);
          success *= (system.evaluateResidual() == 0);
          success *= (system.size() == 27);
          for (IdxT i = 0; i < system.getResidual().getSize(); ++i)
          {
            success *= isEqual(
                system.getResidual().getData()[static_cast<size_t>(i)],
                ZERO<RealT>,
                kTol);
          }
        }

        {
          std::istringstream input(signalSourcedSystemJson());

          auto data  = PhasorDynamics::parseSystemModelData(input);
          success   *= (data.reeca.size() == 1);
          success   *= (data.constant_source.size() == 1);
          success   *= isEqual(
              std::get<RealT>(data.constant_source[0].parameters.at(SourceParams::Sr)),
              static_cast<RealT>(0.6),
              kTol);
          success *= isEqual(
              std::get<RealT>(data.constant_source[0].parameters.at(SourceParams::Si)),
              static_cast<RealT>(0.2),
              kTol);
          success *= (data.constant_source[0].signal_outputs.at(SourceOutputs::sr)
                      == static_cast<IdxT>(10));
          success *= (data.constant_source[0].signal_outputs.at(SourceOutputs::si)
                      == static_cast<IdxT>(11));
          success *= (data.reeca[0].buses.at(Buses::bus) == static_cast<IdxT>(1));
          success *= (data.reeca[0].signal_inputs.at(Inputs::pe) == static_cast<IdxT>(10));
          success *= (data.reeca[0].signal_inputs.at(Inputs::qgen) == static_cast<IdxT>(11));
          success *= (data.reeca[0].signal_outputs.at(Outputs::iqcmd) == static_cast<IdxT>(12));
          success *= (data.reeca[0].signal_outputs.at(Outputs::ipcmd) == static_cast<IdxT>(13));

          PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
          success *= (system.allocate() == 0);
          system.getSignal(12)->init(static_cast<ScalarT>(0.2));
          system.getSignal(13)->init(static_cast<ScalarT>(0.6));
          success *= (system.initialize() == 0);
          success *= isEqual(system.getSignal(10)->read(), static_cast<ScalarT>(0.6), kTol);
          success *= isEqual(system.getSignal(11)->read(), static_cast<ScalarT>(0.2), kTol);
          success *= (system.evaluateResidual() == 0);
          for (IdxT i = 0; i < system.getResidual().getSize(); ++i)
          {
            success *= isEqual(
                system.getResidual().getData()[static_cast<size_t>(i)],
                ZERO<RealT>,
                kTol);
          }
        }

        return success.report(__func__);
      }

    private:
      using Params        = PhasorDynamics::Converter::ReecaParameters;
      using Buses         = PhasorDynamics::Converter::ReecaBuses;
      using Inputs        = PhasorDynamics::Converter::ReecaSignalInputs;
      using Outputs       = PhasorDynamics::Converter::ReecaSignalOutputs;
      using Mon           = PhasorDynamics::Converter::ReecaMonitorableVariables;
      using SourceParams  = PhasorDynamics::ConstantSignalSourceParameters;
      using SourceOutputs = PhasorDynamics::ConstantSignalSourceSignalOutputs;

      static size_t index(PhasorDynamics::Converter::ReecaInternalVariables variable)
      {
        return static_cast<size_t>(variable);
      }

      auto makeMinimalReecaData() -> PhasorDynamics::Converter::ReecaData<RealT, IdxT>
      {
        PhasorDynamics::Converter::ReecaData<RealT, IdxT> data;
        data.device_class            = "Reeca";
        data.disambiguation_string   = "reeca_test";
        data.parameters[Params::mva] = static_cast<RealT>(100.0);
        return data;
      }

      auto makeReecaData() -> PhasorDynamics::Converter::ReecaData<RealT, IdxT>
      {
        auto data = makeMinimalReecaData();
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);

        data.parameters[Params::PfFlag] = static_cast<IdxT>(0);
        data.parameters[Params::VFlag]  = static_cast<IdxT>(1);
        data.parameters[Params::QFlag]  = static_cast<IdxT>(0);
        data.parameters[Params::PFlag]  = static_cast<IdxT>(0);
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
        data.parameters[Params::Iqfrz]  = static_cast<RealT>(0.0);
        data.parameters[Params::Thld]   = static_cast<RealT>(0.0);
        data.parameters[Params::Thld2]  = static_cast<RealT>(0.0);
        data.parameters[Params::Qmax]   = static_cast<RealT>(1.0);
        data.parameters[Params::Qmin]   = static_cast<RealT>(-1.0);
        data.parameters[Params::Kqp]    = static_cast<RealT>(1.0);
        data.parameters[Params::Kqi]    = static_cast<RealT>(0.0);
        data.parameters[Params::Vmax]   = static_cast<RealT>(1.2);
        data.parameters[Params::Vmin]   = static_cast<RealT>(0.8);
        data.parameters[Params::Vref1]  = static_cast<RealT>(0.0);
        data.parameters[Params::Kvp]    = static_cast<RealT>(1.0);
        data.parameters[Params::Kvi]    = static_cast<RealT>(0.0);
        data.parameters[Params::Tiq]    = static_cast<RealT>(0.02);
        data.parameters[Params::Tpord]  = static_cast<RealT>(0.02);
        data.parameters[Params::dPmax]  = static_cast<RealT>(1.0);
        data.parameters[Params::dPmin]  = static_cast<RealT>(-1.0);
        data.parameters[Params::Pmax]   = static_cast<RealT>(1.0);
        data.parameters[Params::Pmin]   = static_cast<RealT>(0.0);
        data.parameters[Params::Imax]   = static_cast<RealT>(2.0);
        data.parameters[Params::Vq1]    = static_cast<RealT>(0.2);
        data.parameters[Params::Iq1]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vq2]    = static_cast<RealT>(0.5);
        data.parameters[Params::Iq2]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vq3]    = static_cast<RealT>(0.75);
        data.parameters[Params::Iq3]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vq4]    = static_cast<RealT>(1.0);
        data.parameters[Params::Iq4]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vp1]    = static_cast<RealT>(0.2);
        data.parameters[Params::Ip1]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vp2]    = static_cast<RealT>(0.5);
        data.parameters[Params::Ip2]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vp3]    = static_cast<RealT>(0.75);
        data.parameters[Params::Ip3]    = static_cast<RealT>(2.0);
        data.parameters[Params::Vp4]    = static_cast<RealT>(1.0);
        data.parameters[Params::Ip4]    = static_cast<RealT>(2.0);

        return data;
      }

      void addAllMonitors(PhasorDynamics::Converter::ReecaData<RealT, IdxT>& data)
      {
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);
      }

      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus, Params param, RealT value)
      {
        auto data              = makeReecaData();
        data.parameters[param] = value;
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> model(&bus, data);
        return model.verify() > 0;
      }

      std::string reecaParamsJson()
      {
        return R"json(
      "params": {
        "mva": 100.0, "PfFlag": 0, "VFlag": 1, "QFlag": 0, "PFlag": 0, "Pqflag": 1,
        "Trv": 0.0, "Tp": 0.02, "Vdip": 0.7, "Vup": 1.2,
        "dbd1": -0.01, "dbd2": 0.01, "kqv": 0.0, "Iql1": -1.0, "Iqh1": 1.0,
        "Iqfrz": 0.0, "Thld": 0.0, "Thld2": 0.0,
        "Qmax": 1.0, "Qmin": -1.0, "Kqp": 1.0, "Kqi": 0.0,
        "Vmax": 1.2, "Vmin": 0.8, "Vref1": 0.0, "Kvp": 1.0, "Kvi": 0.0,
        "Tiq": 0.02, "Tpord": 0.02, "dPmax": 1.0, "dPmin": -1.0,
        "Pmax": 1.0, "Pmin": 0.0, "Imax": 2.0,
        "Vq1": 0.2, "Iq1": 2.0, "Vq2": 0.5, "Iq2": 2.0,
        "Vq3": 0.75, "Iq3": 2.0, "Vq4": 1.0, "Iq4": 2.0,
        "Vp1": 0.2, "Ip1": 2.0, "Vp2": 0.5, "Ip2": 2.0,
        "Vp3": 0.75, "Ip3": 2.0, "Vp4": 1.0, "Ip4": 2.0
      }
)json";
      }

      std::string commandOnlySystemJson()
      {
        return R"json(
{
  "header": {
    "format_version": 0.2,
    "format_revision": 0,
    "case_name": "renewable electrical control",
    "case_description": "REECA parser test",
    "case_comments": ""
  },
  "params": {
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "Bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "params": { "kv": 1.0 } }
  ],
  "signals": [
    { "signal_id": 12, "name": "Iqcmd" },
    { "signal_id": 13, "name": "Ipcmd" }
  ],
  "devices": [
    {
      "class": "Reeca",
      "ports": { "bus": 1, "iqcmd": 12, "ipcmd": 13 },
      "id": "REE1",
)json" + reecaParamsJson()
               +
               R"json(
    }
  ]
}
)json";
      }

      std::string signalSourcedSystemJson()
      {
        return R"json(
{
  "header": {
    "format_version": 0.2,
    "format_revision": 0,
    "case_name": "REECA physics",
    "case_description": "REECA parser and system physics smoke test",
    "case_comments": ""
  },
  "params": {
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "Bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "params": { "kv": 1.0 } }
  ],
  "signals": [
    { "signal_id": 10, "name": "pe fixture" },
    { "signal_id": 11, "name": "qgen fixture" },
    { "signal_id": 12, "name": "iqcmd" },
    { "signal_id": 13, "name": "ipcmd" }
  ],
  "devices": [
    {
      "class": "ConstantSignalSource",
      "ports": { "sr": 10, "si": 11 },
      "id": "PQ",
      "params": { "Sr": 0.6, "Si": 0.2 }
    },
    {
      "class": "Reeca",
      "ports": { "bus": 1, "pe": 10, "qgen": 11, "iqcmd": 12, "ipcmd": 13 },
      "id": "RC",
)json" + reecaParamsJson()
               +
               R"json(
    }
  ]
}
)json";
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      auto makeJacobianData() -> PhasorDynamics::Converter::ReecaData<RealT, IdxT>
      {
        auto data                      = makeReecaData();
        data.parameters[Params::PFlag] = static_cast<IdxT>(1);
        // Distinct VDL limits keep the three linseg coefficients in each of the
        // gq/gp residual rows unequal, which the Enzyme auto-sparsity pass
        // requires to keep structurally identical terms separate.
        data.parameters[Params::Iq1]   = static_cast<RealT>(1.20);
        data.parameters[Params::Iq2]   = static_cast<RealT>(1.15);
        data.parameters[Params::Iq3]   = static_cast<RealT>(1.10);
        data.parameters[Params::Iq4]   = static_cast<RealT>(1.05);
        data.parameters[Params::Ip1]   = static_cast<RealT>(1.20);
        data.parameters[Params::Ip2]   = static_cast<RealT>(1.15);
        data.parameters[Params::Ip3]   = static_cast<RealT>(1.10);
        data.parameters[Params::Ip4]   = static_cast<RealT>(1.05);
        return data;
      }

      DependencyTracking::Variable::DependencyMap pruneJacobianTails(
          DependencyTracking::Variable::DependencyMap dependencies,
          RealT                                       tolerance)
      {
        for (auto it = dependencies.begin(); it != dependencies.end();)
        {
          if (std::abs(it->second) <= tolerance)
          {
            it = dependencies.erase(it);
          }
          else
          {
            ++it;
          }
        }
        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;
        using Var    = PhasorDynamics::Converter::ReecaInternalVariables;
        using Ext    = PhasorDynamics::Converter::ReecaExternalVariables;

        auto data = makeJacobianData();

        PhasorDynamics::Bus<DepVar, IdxT>              bus(DepVar{0.96}, DepVar{0.28});
        PhasorDynamics::Converter::Reeca<DepVar, IdxT> reeca(&bus, data);

        PhasorDynamics::SignalNode<DepVar, IdxT> pe_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> qgen_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> omega_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> ipcmd_node;
        DepVar                                   pe_value{0.8};
        DepVar                                   qgen_value{0.2};
        DepVar                                   omega_value{1.01};
        DepVar                                   iqcmd_value{0.2};
        DepVar                                   ipcmd_value{0.8};
        IdxT                                     pe_index    = static_cast<IdxT>(reeca.size() + bus.size());
        IdxT                                     qgen_index  = static_cast<IdxT>(reeca.size() + bus.size() + 1);
        IdxT                                     omega_index = static_cast<IdxT>(reeca.size() + bus.size() + 2);
        IdxT                                     iqcmd_index = static_cast<IdxT>(reeca.size() + bus.size() + 3);
        IdxT                                     ipcmd_index = static_cast<IdxT>(reeca.size() + bus.size() + 4);

        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        omega_node.set(&omega_value, &omega_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        reeca.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reeca.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        bus.allocate();
        reeca.allocate();
        bus.initialize();

        iqcmd_node.init(DepVar{0.2});
        ipcmd_node.init(DepVar{0.8});
        reeca.initialize();

        for (IdxT i = 0; i < reeca.size(); ++i)
        {
          reeca.y().getData()[static_cast<size_t>(i)].setVariableNumber(i);
          reeca.yp().getData()[static_cast<size_t>(i)].setVariableNumber(i);
        }
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.y().getData()[static_cast<size_t>(i)].setVariableNumber(i + reeca.size());
        }
        reeca.y().setDataUpdated();
        reeca.yp().setDataUpdated();
        bus.y().setDataUpdated();
        pe_value.setVariableNumber(pe_index);
        qgen_value.setVariableNumber(qgen_index);
        omega_value.setVariableNumber(omega_index);

        reeca.evaluateResidual();

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(
            static_cast<size_t>(reeca.size()));
        for (IdxT i = 0; i < reeca.size(); ++i)
        {
          dependencies[static_cast<size_t>(i)] =
              reeca.getResidual().getData()[static_cast<size_t>(i)].getDependencies();
        }

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        using Var = PhasorDynamics::Converter::ReecaInternalVariables;
        using Ext = PhasorDynamics::Converter::ReecaExternalVariables;

        auto data = makeJacobianData();

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(0.96, 0.28);
        PhasorDynamics::Converter::Reeca<ScalarT, IdxT> reeca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        ScalarT                                   pe_value{0.8};
        ScalarT                                   qgen_value{0.2};
        ScalarT                                   omega_value{1.01};
        ScalarT                                   iqcmd_value{0.2};
        ScalarT                                   ipcmd_value{0.8};
        IdxT                                      pe_index    = static_cast<IdxT>(reeca.size() + bus.size());
        IdxT                                      qgen_index  = static_cast<IdxT>(reeca.size() + bus.size() + 1);
        IdxT                                      omega_index = static_cast<IdxT>(reeca.size() + bus.size() + 2);
        IdxT                                      iqcmd_index = static_cast<IdxT>(reeca.size() + bus.size() + 3);
        IdxT                                      ipcmd_index = static_cast<IdxT>(reeca.size() + bus.size() + 4);

        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        omega_node.set(&omega_value, &omega_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        reeca.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reeca.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reeca.getSignals().template attachSignalNode<Ext::OMEGA>(&omega_node);
        reeca.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reeca.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        bus.allocate();
        reeca.allocate();
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + reeca.size());
          bus.setResidualIndex(i, i + reeca.size());
        }

        bus.initialize();
        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.8));
        reeca.initialize();
        reeca.updateTime(0.0, 1.0);

        reeca.evaluateResidual();
        reeca.evaluateJacobian();
        reeca.constructCsr();

        auto model_jacobian = reeca.getCsrJacobian();
        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
