#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ConverterRepcaTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-8);

      TestOutcome constructionAndValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeRepcaData());
        success *= (repca.size()
                    == static_cast<IdxT>(PhasorDynamics::Converter::RepcaInternalVariables::MAXIMUM));
        success *= (repca.getMonitor() != nullptr);
        BranchSignalSet branch_signals;
        branch_signals.attachTo(repca);
        success *= (repca.verify() == 0);

        auto bad_repca                                                        = makeRepcaData();
        bad_repca.parameters[PhasorDynamics::Converter::RepcaParameters::Tfv] = static_cast<RealT>(0.0);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> bad_repca_model(&bus, bad_repca);
        branch_signals.attachTo(bad_repca_model);
        success *= (bad_repca_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome repcaSignalsInitializationAndResidual()
      {
        using Var = PhasorDynamics::Converter::RepcaInternalVariables;
        using Ext = PhasorDynamics::Converter::RepcaExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        ScalarT ibranchr_value{0.6};
        ScalarT ibranchi_value{-0.2};
        ScalarT pbranch_value{0.6};
        ScalarT qbranch_value{0.2};
        ScalarT qext_value{0.0};
        ScalarT pext_value{0.0};
        IdxT    ibranchr_index = 20;
        IdxT    ibranchi_index = 21;
        IdxT    pbranch_index  = 22;
        IdxT    qbranch_index  = 23;
        IdxT    qext_index     = 24;
        IdxT    pext_index     = 25;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ibranchr_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ibranchi_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pbranch_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qbranch_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qext_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pext_node;
        ibranchr_node.set(&ibranchr_value, &ibranchr_index);
        ibranchi_node.set(&ibranchi_value, &ibranchi_index);
        pbranch_node.set(&pbranch_value, &pbranch_index);
        qbranch_node.set(&qbranch_value, &qbranch_index);
        qext_node.set(&qext_value, &qext_index);
        pext_node.set(&pext_value, &pext_index);

        auto                                            data = makeRepcaData();
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ibranchr_node);
        repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ibranchi_node);
        repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&pbranch_node);
        repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&qbranch_node);
        repca.getSignals().template assignSignalNode<Var::QEXT>(&qext_node);
        repca.getSignals().template assignSignalNode<Var::PEXT>(&pext_node);

        success *= (repca.allocate() == 0);
        qext_node.init(static_cast<ScalarT>(0.2));
        pext_node.init(static_cast<ScalarT>(0.6));
        success *= (repca.verify() == 0);
        success *= (repca.initialize() == 0);
        success *= (repca.tagDifferentiable() == 0);
        success *= (repca.evaluateResidual() == 0);

        success *= isEqual(repca.y()[index(Var::PMEAS)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(repca.y()[index(Var::QMEAS)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(repca.y()[index(Var::QEXT)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(repca.y()[index(Var::PEXT)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(qext_node.read(), repca.y()[index(Var::QEXT)], kTol);
        success *= isEqual(pext_node.read(), repca.y()[index(Var::PEXT)], kTol);
        success *= (repca.tag()[index(Var::VMEAS)] == true);
        success *= (repca.tag()[index(Var::PREF)] == false);

        for (size_t i = 0; i < repca.getResidual().size(); ++i)
        {
          success *= isEqual(repca.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(repca.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome rejectsHalfConnectedBranchSignals()
      {
        using Ext = PhasorDynamics::Converter::RepcaExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        ScalarT signal_value{0.6};
        IdxT    signal_index = 30;

        PhasorDynamics::SignalNode<ScalarT, IdxT> signal_node;
        signal_node.set(&signal_value, &signal_index);

        {
          PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeRepcaData());
          repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&signal_node);
          success *= (repca.verify() > 0);
        }

        {
          PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeRepcaData());
          repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&signal_node);
          success *= (repca.verify() > 0);
        }

        {
          PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeRepcaData());
          repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&signal_node);
          success *= (repca.verify() > 0);
        }

        {
          PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeRepcaData());
          repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&signal_node);
          success *= (repca.verify() > 0);
        }

        return success.report(__func__);
      }

      TestOutcome rejectsLineDropCompensationWithoutCurrentSignals()
      {
        using Params = PhasorDynamics::Converter::RepcaParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        auto data                          = makeRepcaData();
        data.parameters[Params::VcompFlag] = static_cast<IdxT>(1);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        BranchSignalSet                                 branch_signals;
        branch_signals.attachPowerTo(repca);

        success *= (repca.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome convertsSystemBaseBranchPowerToComponentBase()
      {
        using Var    = PhasorDynamics::Converter::RepcaInternalVariables;
        using Params = PhasorDynamics::Converter::RepcaParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        BranchSignalSet branch_signals(static_cast<ScalarT>(0.3),
                                       static_cast<ScalarT>(0.1),
                                       static_cast<ScalarT>(0.3),
                                       static_cast<ScalarT>(-0.1),
                                       31);

        auto data                    = makeRepcaData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);
        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        repca.setSystemBase(static_cast<RealT>(60.0), static_cast<RealT>(100.0e6));
        branch_signals.attachTo(repca);

        success *= (repca.allocate() == 0);
        success *= (repca.verify() == 0);
        success *= (repca.initialize() == 0);
        success *= isEqual(repca.y()[index(Var::PMEAS)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(repca.y()[index(Var::QMEAS)], static_cast<ScalarT>(0.2), kTol);
        success *= isEqual(repca.y()[index(Var::PEXT)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(repca.y()[index(Var::QEXT)], static_cast<ScalarT>(0.2), kTol);

        return success.report(__func__);
      }

      TestOutcome mvaZeroUsesSystemBase()
      {
        using Var    = PhasorDynamics::Converter::RepcaInternalVariables;
        using Params = PhasorDynamics::Converter::RepcaParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        BranchSignalSet branch_signals(static_cast<ScalarT>(0.3),
                                       static_cast<ScalarT>(0.1),
                                       static_cast<ScalarT>(0.3),
                                       static_cast<ScalarT>(-0.1),
                                       33);

        auto data                    = makeRepcaData();
        data.parameters[Params::mva] = static_cast<RealT>(0.0);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        repca.setSystemBase(static_cast<RealT>(60.0), static_cast<RealT>(100.0e6));
        branch_signals.attachTo(repca);

        success *= (repca.allocate() == 0);
        success *= (repca.verify() == 0);
        success *= (repca.initialize() == 0);
        success *= isEqual(repca.y()[index(Var::PMEAS)], static_cast<ScalarT>(0.3), kTol);
        success *= isEqual(repca.y()[index(Var::QMEAS)], static_cast<ScalarT>(0.1), kTol);
        success *= isEqual(repca.y()[index(Var::PEXT)], static_cast<ScalarT>(0.3), kTol);
        success *= isEqual(repca.y()[index(Var::QEXT)], static_cast<ScalarT>(0.1), kTol);

        return success.report(__func__);
      }

      TestOutcome rejectsNonSteadyConnectedReferences()
      {
        using Params = PhasorDynamics::Converter::RepcaParameters;
        using Ext    = PhasorDynamics::Converter::RepcaExternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        {
          auto data                   = makeRepcaData();
          data.parameters[Params::Ki] = static_cast<RealT>(1.0);

          ScalarT qref_value{0.5};
          IdxT    qref_index = 35;

          PhasorDynamics::SignalNode<ScalarT, IdxT> qref_node;
          qref_node.set(&qref_value, &qref_index);

          PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
          BranchSignalSet                                 branch_signals;
          branch_signals.attachTo(repca);
          repca.getSignals().template attachSignalNode<Ext::QREF>(&qref_node);

          success *= (repca.allocate() == 0);
          success *= (repca.verify() == 0);
          success *= (repca.initialize() > 0);
        }

        {
          auto data                    = makeRepcaData();
          data.parameters[Params::Kig] = static_cast<RealT>(1.0);

          ScalarT pplantref_value{0.5};
          IdxT    pplantref_index = 36;

          PhasorDynamics::SignalNode<ScalarT, IdxT> pplantref_node;
          pplantref_node.set(&pplantref_value, &pplantref_index);

          PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
          BranchSignalSet                                 branch_signals;
          branch_signals.attachTo(repca);
          repca.getSignals().template attachSignalNode<Ext::PPLANTREF>(&pplantref_node);

          success *= (repca.allocate() == 0);
          success *= (repca.verify() == 0);
          success *= (repca.initialize() > 0);
        }

        return success.report(__func__);
      }

      TestOutcome initializationUsesBranchSignals()
      {
        using Var = PhasorDynamics::Converter::RepcaInternalVariables;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        BranchSignalSet branch_signals(static_cast<ScalarT>(0.1),
                                       static_cast<ScalarT>(-0.4),
                                       static_cast<ScalarT>(0.1),
                                       static_cast<ScalarT>(0.4),
                                       37);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, makeRepcaData());
        branch_signals.attachTo(repca);

        success *= (repca.allocate() == 0);
        success *= (repca.verify() == 0);
        success *= (repca.initialize() == 0);

        success *= isEqual(repca.y()[index(Var::PMEAS)], static_cast<ScalarT>(0.1), kTol);
        success *= isEqual(repca.y()[index(Var::QMEAS)], static_cast<ScalarT>(-0.4), kTol);
        success *= isEqual(repca.y()[index(Var::PEXT)], static_cast<ScalarT>(0.1), kTol);
        success *= isEqual(repca.y()[index(Var::QEXT)], static_cast<ScalarT>(-0.4), kTol);

        success *= (repca.evaluateResidual() == 0);

        for (size_t i = 0; i < repca.getResidual().size(); ++i)
        {
          success *= isEqual(repca.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(repca.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome omittedFrequencyPortsDefaultToZero()
      {
        using Var    = PhasorDynamics::Converter::RepcaInternalVariables;
        using Params = PhasorDynamics::Converter::RepcaParameters;

        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto data                    = makeRepcaData();
        data.parameters[Params::Ddn] = static_cast<RealT>(5.0);
        data.parameters[Params::Dup] = static_cast<RealT>(5.0);
        data.parameters[Params::Kig] = static_cast<RealT>(1.0);

        PhasorDynamics::Converter::Repca<ScalarT, IdxT> repca(&bus, data);
        BranchSignalSet                                 branch_signals;
        branch_signals.attachTo(repca);

        success *= (repca.allocate() == 0);
        success *= (repca.verify() == 0);
        success *= (repca.initialize() == 0);
        success *= (repca.evaluateResidual() == 0);

        const ScalarT pfreq = static_cast<ScalarT>(5.0) * Math::ramp(repca.y()[index(Var::EF)])
                              - static_cast<ScalarT>(5.0) * Math::ramp(-repca.y()[index(Var::EF)]);

        success *= isEqual(repca.y()[index(Var::EF)], static_cast<ScalarT>(0.0), kTol);
        success *= isEqual(pfreq, static_cast<ScalarT>(0.0), kTol);

        for (size_t i = 0; i < repca.getResidual().size(); ++i)
        {
          success *= isEqual(repca.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(repca.yp()[i], static_cast<ScalarT>(0.0), kTol);
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
    "case_name": "renewable plant control",
    "case_description": "REPCA parser test",
    "case_comments": "",
    "freq_base": 60.0,
    "va_base": 100000000.0
  },
  "buses": [
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
  ],
  "signals": [
    { "signal_id": 10, "name": "Qext" },
    { "signal_id": 11, "name": "Pext" },
    { "signal_id": 12, "name": "Ibranchr" },
    { "signal_id": 13, "name": "Ibranchi" },
    { "signal_id": 14, "name": "Pbranch" },
    { "signal_id": 15, "name": "Qbranch" }
  ],
  "devices": [
    {
      "class": "Regca",
      "ports": { "bus": 1, "ibranchr": 12, "ibranchi": 13, "pbranch": 14, "qbranch": 15 },
      "id": "RG1",
      "params": {
        "P0": 0.0,
        "Q0": 0.0,
        "mva": 100.0,
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
      }
    },
    {
      "class": "Repca",
      "ports": {
        "bus": 1,
        "ibranchr": 12,
        "ibranchi": 13,
        "pbranch": 14,
        "qbranch": 15,
        "qext": 10,
        "pext": 11
      },
      "id": "REP1",
      "params": {
        "mva": 100.0, "VcompFlag": 0, "RefFlag": 0, "Freqflag": 1,
        "Tfltr": 0.02, "Tft": 0.0, "Tfv": 0.02, "Tp": 0.02, "Tlag": 0.0,
        "Vfrz": 0.7, "Rc": 0.0, "Xc": 0.0, "Kc": 0.0,
        "dbdlow": -0.01, "dbdupper": 0.01, "emax": 0.1, "emin": -0.1,
        "Kp": 1.0, "Ki": 0.0, "Qmax": 1.0, "Qmin": -1.0,
        "fdbd1": -0.01, "fdbd2": 0.01, "Ddn": 0.0, "Dup": 0.0,
        "femax": 0.1, "femin": -0.1, "Kpg": 1.0, "Kig": 0.0,
        "Pmax": 1.0, "Pmin": 0.0
      }
    }
  ]
}
)json");

        auto data  = PhasorDynamics::parseSystemModelData(input);
        success   *= (data.regca.size() == 1);
        success   *= (data.repca.size() == 1);
        success   *= (data.repca[0].ports.at(PhasorDynamics::Converter::RepcaPorts::pbranch) == 14);
        success   *= (data.repca[0].ports.at(PhasorDynamics::Converter::RepcaPorts::qbranch) == 15);
        success   *= (std::get<RealT>(
                        data.repca[0].parameters.at(PhasorDynamics::Converter::RepcaParameters::mva))
                    == static_cast<RealT>(100.0));

        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() > 0);

        return success.report(__func__);
      }

    private:
      static size_t index(PhasorDynamics::Converter::RepcaInternalVariables variable)
      {
        return static_cast<size_t>(variable);
      }

      struct BranchSignalSet
      {
        ScalarT ibranchr_value;
        ScalarT ibranchi_value;
        ScalarT pbranch_value;
        ScalarT qbranch_value;
        IdxT    ibranchr_index;
        IdxT    ibranchi_index;
        IdxT    pbranch_index;
        IdxT    qbranch_index;

        PhasorDynamics::SignalNode<ScalarT, IdxT> ibranchr_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ibranchi_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pbranch_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qbranch_node;

        BranchSignalSet(ScalarT pbranch     = static_cast<ScalarT>(0.6),
                        ScalarT qbranch     = static_cast<ScalarT>(0.2),
                        ScalarT ibranchr    = static_cast<ScalarT>(0.6),
                        ScalarT ibranchi    = static_cast<ScalarT>(-0.2),
                        IdxT    first_index = 20)
          : ibranchr_value(ibranchr),
            ibranchi_value(ibranchi),
            pbranch_value(pbranch),
            qbranch_value(qbranch),
            ibranchr_index(first_index),
            ibranchi_index(first_index + 1),
            pbranch_index(first_index + 2),
            qbranch_index(first_index + 3)
        {
          ibranchr_node.set(&ibranchr_value, &ibranchr_index);
          ibranchi_node.set(&ibranchi_value, &ibranchi_index);
          pbranch_node.set(&pbranch_value, &pbranch_index);
          qbranch_node.set(&qbranch_value, &qbranch_index);
        }

        template <class RepcaT>
        void attachTo(RepcaT& repca)
        {
          using Ext = PhasorDynamics::Converter::RepcaExternalVariables;

          repca.getSignals().template attachSignalNode<Ext::IBRANCHR>(&ibranchr_node);
          repca.getSignals().template attachSignalNode<Ext::IBRANCHI>(&ibranchi_node);
          repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&pbranch_node);
          repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&qbranch_node);
        }

        template <class RepcaT>
        void attachPowerTo(RepcaT& repca)
        {
          using Ext = PhasorDynamics::Converter::RepcaExternalVariables;

          repca.getSignals().template attachSignalNode<Ext::PBRANCH>(&pbranch_node);
          repca.getSignals().template attachSignalNode<Ext::QBRANCH>(&qbranch_node);
        }
      };

      auto makeRepcaData() -> PhasorDynamics::Converter::RepcaData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Converter::RepcaParameters;
        using Mon    = PhasorDynamics::Converter::RepcaMonitorableVariables;

        PhasorDynamics::Converter::RepcaData<RealT, IdxT> data;
        data.device_class          = "Repca";
        data.disambiguation_string = "repca_test";
        data.monitored_variables.insert(Mon::qext);
        data.monitored_variables.insert(Mon::pext);

        data.parameters[Params::mva]       = static_cast<RealT>(100.0);
        data.parameters[Params::VcompFlag] = static_cast<IdxT>(0);
        data.parameters[Params::RefFlag]   = static_cast<IdxT>(0);
        data.parameters[Params::Freqflag]  = static_cast<IdxT>(1);
        data.parameters[Params::Tfltr]     = static_cast<RealT>(0.02);
        data.parameters[Params::Tft]       = static_cast<RealT>(0.0);
        data.parameters[Params::Tfv]       = static_cast<RealT>(0.02);
        data.parameters[Params::Tp]        = static_cast<RealT>(0.02);
        data.parameters[Params::Tlag]      = static_cast<RealT>(0.0);
        data.parameters[Params::Vfrz]      = static_cast<RealT>(0.7);
        data.parameters[Params::Rc]        = static_cast<RealT>(0.0);
        data.parameters[Params::Xc]        = static_cast<RealT>(0.0);
        data.parameters[Params::Kc]        = static_cast<RealT>(0.0);
        data.parameters[Params::dbdlow]    = static_cast<RealT>(-0.01);
        data.parameters[Params::dbdupper]  = static_cast<RealT>(0.01);
        data.parameters[Params::emax]      = static_cast<RealT>(0.1);
        data.parameters[Params::emin]      = static_cast<RealT>(-0.1);
        data.parameters[Params::Kp]        = static_cast<RealT>(1.0);
        data.parameters[Params::Ki]        = static_cast<RealT>(0.0);
        data.parameters[Params::Qmax]      = static_cast<RealT>(1.0);
        data.parameters[Params::Qmin]      = static_cast<RealT>(-1.0);
        data.parameters[Params::fdbd1]     = static_cast<RealT>(-0.01);
        data.parameters[Params::fdbd2]     = static_cast<RealT>(0.01);
        data.parameters[Params::Ddn]       = static_cast<RealT>(0.0);
        data.parameters[Params::Dup]       = static_cast<RealT>(0.0);
        data.parameters[Params::femax]     = static_cast<RealT>(0.1);
        data.parameters[Params::femin]     = static_cast<RealT>(-0.1);
        data.parameters[Params::Kpg]       = static_cast<RealT>(1.0);
        data.parameters[Params::Kig]       = static_cast<RealT>(0.0);
        data.parameters[Params::Pmax]      = static_cast<RealT>(1.0);
        data.parameters[Params::Pmin]      = static_cast<RealT>(0.0);

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
