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

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, makeReecbData());
        success *=
            (reecb.size()
             == static_cast<IdxT>(PhasorDynamics::Converter::ReecbInternalVariables::MAXIMUM));
        success *= (reecb.getMonitor() != nullptr);
        success *= (reecb.verify() == 0);

        auto wide_band_reecb = makeReecbData();
        wide_band_reecb.parameters[PhasorDynamics::Converter::ReecbParameters::Vdip] =
            static_cast<RealT>(-99.0);
        wide_band_reecb.parameters[PhasorDynamics::Converter::ReecbParameters::Vup] =
            static_cast<RealT>(99.0);
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> wide_band_reecb_model(
            &bus,
            wide_band_reecb);
        success *= (wide_band_reecb_model.verify() == 0);

        auto bad_reecb_band = makeReecbData();
        bad_reecb_band.parameters[PhasorDynamics::Converter::ReecbParameters::Vdip] =
            static_cast<RealT>(1.2);
        bad_reecb_band.parameters[PhasorDynamics::Converter::ReecbParameters::Vup] =
            static_cast<RealT>(1.2);
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> bad_reecb_band_model(
            &bus,
            bad_reecb_band);
        success *= (bad_reecb_band_model.verify() > 0);

        auto bad_reecb = makeReecbData();
        bad_reecb.parameters[PhasorDynamics::Converter::ReecbParameters::Imax] =
            static_cast<RealT>(-1.0);
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> bad_reecb_model(&bus, bad_reecb);
        success *= (bad_reecb_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome reecbSignalsInitializationAndResidual()
      {
        using Var = PhasorDynamics::Converter::ReecbInternalVariables;
        using Ext = PhasorDynamics::Converter::ReecbExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ScalarT pe_value{0.6};
        ScalarT qgen_value{0.2};
        ScalarT iqcmd_value{0.2};
        ScalarT ipcmd_value{0.6};
        IdxT    pe_index    = 20;
        IdxT    qgen_index  = 21;
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        auto data                                                          = makeReecbData();
        data.parameters[PhasorDynamics::Converter::ReecbParameters::QFlag] = static_cast<IdxT>(1);
        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, data);
        reecb.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reecb.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        success *= (reecb.verify() == 0);
        iqcmd_node.init(qgen_value);
        ipcmd_node.init(pe_value);
        success *= (reecb.initialize() == 0);
        success *= (reecb.tagDifferentiable() == 0);
        success *= (reecb.evaluateResidual() == 0);

        success *= isEqual(reecb.y()[index(Var::VMEAS)], static_cast<ScalarT>(1.0), kTol);
        success *= isEqual(reecb.y()[index(Var::PMEAS)], pe_value, kTol);
        success *= isEqual(reecb.y()[index(Var::QREF)], qgen_value, kTol);
        success *= isEqual(iqcmd_node.read(), reecb.y()[index(Var::IQCMD)], kTol);
        success *= isEqual(ipcmd_node.read(), reecb.y()[index(Var::IPCMD)], kTol);
        success *= (reecb.tag()[index(Var::VMEAS)] == false);
        success *= (reecb.tag()[index(Var::PMEAS)] == true);

        for (size_t i = 0; i < reecb.getResidual().size(); ++i)
        {
          success *= isEqual(reecb.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reecb.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome rejectsHalfConnectedElectricalFeedback()
      {
        using Ext = PhasorDynamics::Converter::ReecbExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        ScalarT signal_value{0.6};
        IdxT    signal_index = 24;

        PhasorDynamics::SignalNode<ScalarT, IdxT> signal_node;
        signal_node.set(&signal_value, &signal_index);

        {
          PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, makeReecbData());
          reecb.getSignals().template attachSignalNode<Ext::PE>(&signal_node);
          success *= (reecb.verify() > 0);
        }

        {
          PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, makeReecbData());
          reecb.getSignals().template attachSignalNode<Ext::QGEN>(&signal_node);
          success *= (reecb.verify() > 0);
        }

        return success.report(__func__);
      }

      TestOutcome reecbCommandSignalInitialization()
      {
        using Var = PhasorDynamics::Converter::ReecbInternalVariables;

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

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, makeReecbData());
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        success *= (reecb.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.2));
        ipcmd_node.init(static_cast<ScalarT>(0.6));
        success *= (reecb.initialize() == 0);
        success *= (reecb.tagDifferentiable() == 0);
        success *= (reecb.evaluateResidual() == 0);

        success *= isEqual(reecb.y()[index(Var::PMEAS)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reecb.y()[index(Var::QREF)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(reecb.y()[index(Var::PORD)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reecb.y()[index(Var::IPCMD)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(reecb.y()[index(Var::IQCMD)], static_cast<ScalarT>(0.2), kTol);

        for (size_t i = 0; i < reecb.getResidual().size(); ++i)
        {
          success *= isEqual(reecb.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reecb.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome reecbElectricalFeedbackUsesMvaBase()
      {
        using Var    = PhasorDynamics::Converter::ReecbInternalVariables;
        using Params = PhasorDynamics::Converter::ReecbParameters;
        using Ext    = PhasorDynamics::Converter::ReecbExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeReecbData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        ScalarT pe_value{0.25};
        ScalarT qgen_value{0.05};
        ScalarT iqcmd_value{0.1};
        ScalarT ipcmd_value{0.5};
        IdxT    pe_index    = 25;
        IdxT    qgen_index  = 26;
        IdxT    iqcmd_index = 27;
        IdxT    ipcmd_index = 28;

        PhasorDynamics::SignalNode<ScalarT, IdxT> pe_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qgen_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        pe_node.set(&pe_value, &pe_index);
        qgen_node.set(&qgen_value, &qgen_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, data);
        reecb.getSignals().template attachSignalNode<Ext::PE>(&pe_node);
        reecb.getSignals().template attachSignalNode<Ext::QGEN>(&qgen_node);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        success *= (reecb.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.1));
        ipcmd_node.init(static_cast<ScalarT>(0.5));
        success *= (reecb.initialize() == 0);
        success *= (reecb.evaluateResidual() == 0);

        success *= isEqual(reecb.y()[index(Var::PMEAS)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(reecb.y()[index(Var::QREF)], static_cast<ScalarT>(0.1), kTol);
        success *= isEqual(reecb.y()[index(Var::PORD)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(reecb.y()[index(Var::IPCMD)], static_cast<ScalarT>(0.5), kTol);
        success *= isEqual(reecb.y()[index(Var::IQCMD)], static_cast<ScalarT>(0.1), kTol);

        return success.report(__func__);
      }

      TestOutcome reecbReferenceFallbackAtAngle()
      {
        using Var    = PhasorDynamics::Converter::ReecbInternalVariables;
        using Params = PhasorDynamics::Converter::ReecbParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        auto data                       = makeReecbData();
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

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, data);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        success *= (reecb.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(-0.2));
        ipcmd_node.init(static_cast<ScalarT>(1.6));
        success *= (reecb.initialize() == 0);
        success *= (reecb.tagDifferentiable() == 0);
        success *= (reecb.evaluateResidual() == 0);

        success *= isEqual(reecb.y()[index(Var::PORD)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(reecb.y()[index(Var::QREF)], static_cast<ScalarT>(-0.2), kTol);
        success *= isEqual(reecb.y()[index(Var::IPCMD)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(reecb.y()[index(Var::IQCMD)], static_cast<ScalarT>(-0.2), kTol);

        for (size_t i = 0; i < reecb.getResidual().size(); ++i)
        {
          success *= isEqual(reecb.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(reecb.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome initializesSaturatedStartConsistently()
      {
        using Var    = PhasorDynamics::Converter::ReecbInternalVariables;
        using Params = PhasorDynamics::Converter::ReecbParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.13, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                       = makeReecbData();
        data.parameters[Params::QFlag]  = static_cast<IdxT>(0);
        data.parameters[Params::Pqflag] = static_cast<IdxT>(0);
        data.parameters[Params::Vmin]   = static_cast<RealT>(0.9);
        data.parameters[Params::Vmax]   = static_cast<RealT>(1.05);
        data.parameters[Params::Kvp]    = static_cast<RealT>(10.0);
        data.parameters[Params::Kvi]    = static_cast<RealT>(60.0);
        data.parameters[Params::Vup]    = static_cast<RealT>(99.0);
        data.parameters[Params::Vdip]   = static_cast<RealT>(-99.0);
        data.parameters[Params::Imax]   = static_cast<RealT>(1.1);

        ScalarT iqcmd_value{0.15 / 1.13};
        ScalarT ipcmd_value{0.5 / 1.13};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, data);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        success *= (reecb.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.15 / 1.13));
        ipcmd_node.init(static_cast<ScalarT>(0.5 / 1.13));
        success *= (reecb.initialize() == 0);

        const ScalarT piv_arg = static_cast<ScalarT>(10.0) * reecb.y()[index(Var::EPIV)]
                                + reecb.y()[index(Var::XPIV)];
        success *= (piv_arg < -reecb.y()[index(Var::IQMAX)]);

        success *= (reecb.evaluateResidual() == 0);
        for (size_t i = 0; i < reecb.getResidual().size(); ++i)
        {
          success *= isEqual(reecb.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome reecbSaturatedVoltagePiJacobianFinite()
      {
        using Var    = PhasorDynamics::Converter::ReecbInternalVariables;
        using Params = PhasorDynamics::Converter::ReecbParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.13, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                       = makeReecbData();
        data.parameters[Params::QFlag]  = static_cast<IdxT>(0);
        data.parameters[Params::Pqflag] = static_cast<IdxT>(0);
        data.parameters[Params::Vmin]   = static_cast<RealT>(0.9);
        data.parameters[Params::Vmax]   = static_cast<RealT>(1.05);
        data.parameters[Params::Kvp]    = static_cast<RealT>(10.0);
        data.parameters[Params::Kvi]    = static_cast<RealT>(60.0);
        data.parameters[Params::Vup]    = static_cast<RealT>(99.0);
        data.parameters[Params::Vdip]   = static_cast<RealT>(-99.0);
        data.parameters[Params::Imax]   = static_cast<RealT>(1.1);

        ScalarT iqcmd_value{0.15 / 1.13};
        ScalarT ipcmd_value{0.5 / 1.13};
        IdxT    iqcmd_index = 22;
        IdxT    ipcmd_index = 23;

        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);

        PhasorDynamics::Converter::Reecb<ScalarT, IdxT> reecb(&bus, data);
        reecb.getSignals().template assignSignalNode<Var::IQCMD>(&iqcmd_node);
        reecb.getSignals().template assignSignalNode<Var::IPCMD>(&ipcmd_node);

        success *= (reecb.allocate() == 0);
        success *= (reecb.verify() == 0);
        iqcmd_node.init(static_cast<ScalarT>(0.15 / 1.13));
        ipcmd_node.init(static_cast<ScalarT>(0.5 / 1.13));
        success *= (reecb.initialize() == 0);
        success *= (reecb.tagDifferentiable() == 0);
        success *= (reecb.evaluateResidual() == 0);
        success *= (reecb.evaluateJacobian() == 0);

        const ScalarT piv_arg = static_cast<ScalarT>(10.0) * reecb.y()[index(Var::EPIV)]
                                + reecb.y()[index(Var::XPIV)];
        success *= (piv_arg < -reecb.y()[index(Var::IQMAX)]);

        auto        jacobian_entries  = reecb.getJacobian().getEntries(false);
        const auto& values            = std::get<2>(jacobian_entries);
        success                      *= (!values.empty());
        for (const auto value : values)
        {
          success *= std::isfinite(value);
        }

        return success.report(__func__);
      }
#endif

      TestOutcome jsonParseAndSystemAssembly()
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
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
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
        success   *= (std::get<IdxT>(
                        data.reecb[0].parameters.at(PhasorDynamics::Converter::ReecbParameters::Pqflag))
                    == static_cast<IdxT>(1));
        success   *= (std::get<RealT>(
                        data.reecb[0].parameters.at(PhasorDynamics::Converter::ReecbParameters::mva))
                    == static_cast<RealT>(100.0));

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() == 27);

        return success.report(__func__);
      }

    private:
      static size_t index(PhasorDynamics::Converter::ReecbInternalVariables variable)
      {
        return static_cast<size_t>(variable);
      }

      auto makeReecbData() -> PhasorDynamics::Converter::ReecbData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Converter::ReecbParameters;
        using Mon    = PhasorDynamics::Converter::ReecbMonitorableVariables;

        PhasorDynamics::Converter::ReecbData<RealT, IdxT> data;
        data.device_class          = "Reecb";
        data.disambiguation_string = "reecb_test";
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);

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
