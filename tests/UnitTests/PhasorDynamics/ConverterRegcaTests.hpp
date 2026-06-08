#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ConverterRegcaTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterRegcaTests()  = default;
      ~ConverterRegcaTests() = default;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-9);

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
        success *= invalidParameterCase(bus, Params::Tg, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::TM, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Rpmax, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::Rqmin, static_cast<RealT>(0.0));
        success *= invalidParameterCase(bus, Params::IL1, static_cast<RealT>(-0.1));
        success *= invalidParameterCase(bus, Params::VL1, static_cast<RealT>(0.3));
        success *= invalidParameterCase(bus, Params::VA1, static_cast<RealT>(0.3));
        success *= invalidParameterCase(bus, Params::Vhvmax, static_cast<RealT>(0.0));

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
        const ScalarT ipcmd0 = toComponentBase(p0 / vt, mva) / lvacm;
        const ScalarT iqcmd0 = toComponentBase(q0 / vt, mva);
        const ScalarT qnet0  = iqcmd0;
        const ScalarT ir0    = (lvacm * ipcmd0 * vr + qnet0 * vi) / vt;
        const ScalarT ii0    = (lvacm * ipcmd0 * vi - qnet0 * vr) / vt;

        success *= scalarMatches(regca.y()[index(Vars::VM)], vt, "VM");
        success *= scalarMatches(regca.y()[index(Vars::VT)], vt, "VT");
        success *= scalarMatches(regca.y()[index(Vars::IP)], ipcmd0, "IP");
        success *= scalarMatches(regca.y()[index(Vars::IQ)], iqcmd0, "IQ");
        success *= scalarMatches(regca.y()[index(Vars::IR)], ir0, "IR");
        success *= scalarMatches(regca.y()[index(Vars::II)], ii0, "II");
        success *= scalarMatches(regca.y()[index(Vars::PBR)], p0, "PBR");
        success *= scalarMatches(regca.y()[index(Vars::QBR)], q0, "QBR");
        success *= scalarMatches(ipcmd_node.read(), ipcmd0, "ipcmd signal");
        success *= scalarMatches(iqcmd_node.read(), iqcmd0, "iqcmd signal");
        success *= scalarMatches(ir_node.read(), ir0, "ir signal");
        success *= scalarMatches(ii_node.read(), ii0, "ii signal");
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

        success *= (regca.allocate() == 0);
        success *= (regca.verify() == 0);
        success *= (regca.initialize() == 0);
        success *= (regca.evaluateResidual() == 0);

        success *= isEqual(regca.y()[index(Vars::IP)], static_cast<ScalarT>(0.6), kTol);
        success *= isEqual(regca.y()[index(Vars::IQ)], static_cast<ScalarT>(0.2), kTol);
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
        success *= (regca.evaluateResidual() == 0);

        const ScalarT tg       = static_cast<RealT>(0.02);
        const ScalarT ip_error = ipcmd_value - regca.y()[index(Vars::IP)];
        const ScalarT iq_error = iqcmd_value - regca.y()[index(Vars::IQ)];

        const ScalarT expected_ip_residual =
            tg * regca.y()[index(Vars::LP)]
            + Math::ramp(ip_error - tg * regca.y()[index(Vars::LP)])
            - Math::ramp(ip_error - tg * regca.y()[index(Vars::UP)]);
        const ScalarT expected_iq_residual =
            iq_error + Math::ramp(tg * static_cast<RealT>(-999.0) - iq_error);

        success *= isEqual(regca.getResidual()[index(Vars::IP)], expected_ip_residual, kTol);
        success *= isEqual(regca.getResidual()[index(Vars::IQ)], expected_iq_residual, kTol);
        success *= (regca.getResidual()[index(Vars::IP)] > ZERO<RealT>);
        success *= (regca.getResidual()[index(Vars::IQ)] > ZERO<RealT>);

        return success.report(__func__);
      }

      TestOutcome invalidInitialization()
      {
        TestStatus success = true;

        auto data = makeData();

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(1.3, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
          success *= (regca.allocate() == 0);
          success *= (regca.initialize() > 0);
        }

        {
          PhasorDynamics::Bus<ScalarT, IdxT> bus(0.0, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
          success *= (regca.allocate() == 0);
          success *= (regca.initialize() > 0);
        }

        {
          auto low_voltage_data                   = data;
          low_voltage_data.parameters[Params::P0] = static_cast<RealT>(0.5);

          PhasorDynamics::Bus<ScalarT, IdxT> bus(0.2, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, low_voltage_data);
          success *= (regca.allocate() == 0);
          success *= (regca.initialize() > 0);
        }

        {
          auto low_voltage_data                   = data;
          low_voltage_data.parameters[Params::P0] = static_cast<RealT>(0.0);
          low_voltage_data.parameters[Params::Q0] = static_cast<RealT>(0.1);

          PhasorDynamics::Bus<ScalarT, IdxT> bus(0.2, 0.0);
          bus.allocate();
          bus.initialize();

          PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, low_voltage_data);
          success *= (regca.allocate() == 0);
          success *= (regca.initialize() == 0);
        }

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
        success *= (regca.initialize() > 0);

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
        success *= (regca.evaluateResidual() == 0);

        success *= isEqual(bus.Ir(), toSystemBase(regca.y()[index(Vars::IR)], mva), kTol);
        success *= isEqual(bus.Ii(), toSystemBase(regca.y()[index(Vars::II)], mva), kTol);
        success *= isEqual(bus.Ir(), static_cast<ScalarT>(0.26), kTol);
        success *= isEqual(bus.Ii(), static_cast<ScalarT>(0.32), kTol);

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

        setResidualState(regca);
        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);

        const auto expected  = expectedResidualForState(ipcmd_value, iqcmd_value);
        success             *= vectorMatches(regca.getResidual(), expected, "REGCA residual");

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

        regca.y()[index(Vars::VT)] = static_cast<ScalarT>(1.35);
        regca.y()[index(Vars::IQEXTRA)] =
            Math::ramp(regca.y()[index(Vars::VT)] - static_cast<RealT>(1.3));

        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);
        success *= scalarMatches(regca.getResidual()[index(Vars::IQEXTRA)],
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
          success *= (regca.evaluateResidual() == 0);

          const ScalarT iq_error = iqcmd_value - regca.y()[index(Vars::IQ)];
          const ScalarT expected =
              iq_error - Math::ramp(iq_error - static_cast<RealT>(0.2) * static_cast<RealT>(0.5));

          success *= scalarMatches(regca.getResidual()[index(Vars::IQ)],
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
          success *= (regca.evaluateResidual() == 0);

          const ScalarT ip       = regca.y()[index(Vars::IP)];
          const ScalarT sigma_ip = Math::sigmoid(ip);
          const RealT   rpmax    = static_cast<RealT>(0.7);
          const RealT   mp       = static_cast<RealT>(100.0) * rpmax;
          const ScalarT expected = mp * (ONE<RealT> - sigma_ip) + rpmax * sigma_ip;

          success *= scalarMatches(regca.y()[index(Vars::UP)], expected, "UP without LVPL");
          success *= scalarMatches(regca.getResidual()[index(Vars::UP)],
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
      "v_base": 1.0
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
        success                *= (regca_data.ports.at(PhasorDynamics::Converter::RegcaPorts::bus)
                    == 1);
        success                *= (std::get_if<double>(&regca_data.parameters.at(Params::P0))
                    != nullptr);
        success                *= (std::get_if<double>(&regca_data.parameters.at(Params::Q0))
                    != nullptr);
        success                *= (std::get_if<size_t>(&regca_data.parameters.at(Params::mva))
                    != nullptr);
        success                *= (std::get_if<bool>(&regca_data.parameters.at(Params::sL))
                    != nullptr);

        PhasorDynamics::SystemModel<double, size_t> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.tagDifferentiable() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.evaluateJacobian() == 0);
        success *= (system.size() == 14);
        success *= isEqual(system.getResidual()[0], 0.0, static_cast<double>(kTol));
        success *= isEqual(system.getResidual()[1], 0.0, static_cast<double>(kTol));
        for (size_t i = 2; i < system.getResidual().size(); ++i)
        {
          success *= isEqual(system.getResidual()[i], 0.0, static_cast<double>(kTol));
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
        using Ports = PhasorDynamics::Converter::RegcaPorts;
        using Mon   = PhasorDynamics::Converter::RegcaMonitorableVariables;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> data;
        data.device_class          = "Regca";
        data.disambiguation_string = "regca_test";
        data.ports[Ports::bus]     = 1;
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

      ScalarT activeLower(ScalarT ip, RealT rpmax) const
      {
        const RealT mp = static_cast<RealT>(100.0) * rpmax;
        return -rpmax - (mp - rpmax) * Math::sigmoid(ip);
      }

      ScalarT activeUpper(ScalarT ip, ScalarT il, RealT rpmax, bool use_lvpl) const
      {
        const RealT   mp       = static_cast<RealT>(100.0) * rpmax;
        const ScalarT sigma_ip = Math::sigmoid(ip);
        ScalarT       output   = mp * (ONE<RealT> - sigma_ip) + rpmax * sigma_ip;
        if (use_lvpl)
        {
          output = mp * (ONE<RealT> - sigma_ip) + rpmax * sigma_ip * Math::sigmoid(il - ip);
        }
        return output;
      }

      bool allResidualsZero(PhasorDynamics::Converter::Regca<ScalarT, IdxT>& regca) const
      {
        bool success = true;
        for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
        {
          if (!isEqual(regca.getResidual()[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "REGCA residual row " << i << " is " << regca.getResidual()[i] << "\n";
            success = false;
          }
          if (!isEqual(regca.yp()[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "REGCA derivative row " << i << " is " << regca.yp()[i] << "\n";
            success = false;
          }
        }
        return success;
      }

      bool vectorMatches(const std::vector<ScalarT>& actual,
                         const std::vector<ScalarT>& expected,
                         const char*                 label) const
      {
        bool       success = (actual.size() == expected.size());
        const auto n       = std::min(actual.size(), expected.size());
        for (size_t i = 0; i < n; ++i)
        {
          if (!isEqual(actual[i], expected[i], kTol))
          {
            std::cout << label << " mismatch at row " << i << ": "
                      << actual[i] << " != " << expected[i] << "\n";
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

        std::cout << label << " mismatch: " << actual << " != " << expected << "\n";
        return false;
      }

      void setResidualState(PhasorDynamics::Converter::Regca<ScalarT, IdxT>& regca)
      {
        regca.y()[index(Vars::VM)]      = static_cast<ScalarT>(0.86);
        regca.y()[index(Vars::IQ)]      = static_cast<ScalarT>(-0.2);
        regca.y()[index(Vars::IP)]      = static_cast<ScalarT>(0.85);
        regca.y()[index(Vars::VT)]      = static_cast<ScalarT>(0.98);
        regca.y()[index(Vars::II)]      = static_cast<ScalarT>(0.18);
        regca.y()[index(Vars::IQEXTRA)] = static_cast<ScalarT>(0.03);
        regca.y()[index(Vars::IL)]      = static_cast<ScalarT>(0.72);
        regca.y()[index(Vars::IR)]      = static_cast<ScalarT>(0.5);
        regca.y()[index(Vars::LP)]      = static_cast<ScalarT>(-0.4);
        regca.y()[index(Vars::UP)]      = static_cast<ScalarT>(0.3);
        regca.y()[index(Vars::PBR)]     = static_cast<ScalarT>(0.52);
        regca.y()[index(Vars::QBR)]     = static_cast<ScalarT>(-0.046);

        regca.yp()[index(Vars::VM)] = static_cast<ScalarT>(0.01);
        regca.yp()[index(Vars::IQ)] = static_cast<ScalarT>(-0.02);
        regca.yp()[index(Vars::IP)] = static_cast<ScalarT>(0.03);
      }

      std::vector<ScalarT> expectedResidualForState(ScalarT ipcmd, ScalarT iqcmd) const
      {
        constexpr RealT   tg    = static_cast<RealT>(0.2);
        constexpr RealT   tm    = static_cast<RealT>(0.4);
        constexpr RealT   rqmin = static_cast<RealT>(-0.6);
        constexpr RealT   rpmax = static_cast<RealT>(0.7);
        constexpr RealT   vl0   = static_cast<RealT>(0.4);
        constexpr RealT   vl1   = static_cast<RealT>(0.9);
        constexpr RealT   va0   = static_cast<RealT>(0.4);
        constexpr RealT   va1   = static_cast<RealT>(0.9);
        constexpr RealT   il1   = static_cast<RealT>(1.1);
        constexpr ScalarT vr    = static_cast<ScalarT>(0.95);
        constexpr ScalarT vi    = static_cast<ScalarT>(0.25);
        constexpr ScalarT vm    = static_cast<ScalarT>(0.86);
        constexpr ScalarT iq    = static_cast<ScalarT>(-0.2);
        constexpr ScalarT ip    = static_cast<ScalarT>(0.85);
        constexpr ScalarT vt    = static_cast<ScalarT>(0.98);
        constexpr ScalarT ii    = static_cast<ScalarT>(0.18);
        constexpr ScalarT iqext = static_cast<ScalarT>(0.03);
        constexpr ScalarT il    = static_cast<ScalarT>(0.72);
        constexpr ScalarT ir    = static_cast<ScalarT>(0.5);
        constexpr ScalarT lp    = static_cast<ScalarT>(-0.4);
        constexpr ScalarT up    = static_cast<ScalarT>(0.3);
        constexpr ScalarT pbr   = static_cast<ScalarT>(0.52);
        constexpr ScalarT qbr   = static_cast<ScalarT>(-0.046);
        constexpr ScalarT vmdot = static_cast<ScalarT>(0.01);
        constexpr ScalarT iqdot = static_cast<ScalarT>(-0.02);
        constexpr ScalarT ipdot = static_cast<ScalarT>(0.03);

        const ScalarT iq_error = iqcmd - iq;
        const ScalarT ip_error = ipcmd - ip;
        const ScalarT lvacm    = Math::linseg(vt, va0, va1, ONE<RealT>);
        const ScalarT qnet     = iq - iqext;

        std::vector<ScalarT> expected(static_cast<size_t>(Vars::MAXIMUM),
                                      static_cast<ScalarT>(0.0));
        expected[index(Vars::VM)] = -tm * vmdot - vm + vt;
        expected[index(Vars::IQ)] = -tg * iqdot + iq_error + Math::ramp(tg * rqmin - iq_error);
        expected[index(Vars::IP)] =
            -tg * ipdot + tg * lp + Math::ramp(ip_error - tg * lp) - Math::ramp(ip_error - tg * up);
        expected[index(Vars::VT)]      = -vt * vt + vr * vr + vi * vi;
        expected[index(Vars::II)]      = -vt * ii + lvacm * ip * vi - qnet * vr;
        expected[index(Vars::IQEXTRA)] = -iqext + Math::ramp(vt - static_cast<RealT>(1.3));
        expected[index(Vars::IL)]      = -il + Math::linseg(vm, vl0, vl1, il1);
        expected[index(Vars::IR)]      = -vt * ir + lvacm * ip * vr + qnet * vi;
        expected[index(Vars::LP)]      = -lp + activeLower(ip, rpmax);
        expected[index(Vars::UP)]      = -up + activeUpper(ip, il, rpmax, true);
        expected[index(Vars::PBR)]     = -pbr + vr * ir + vi * ii;
        expected[index(Vars::QBR)]     = -qbr + vi * ir - vr * ii;
        return expected;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      void setResidualStateDep(
          PhasorDynamics::Converter::Regca<DependencyTracking::Variable, IdxT>& regca,
          PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>&              bus)
      {
        bus.y()[0].setValue(0.95);
        bus.y()[1].setValue(0.25);

        regca.y()[index(Vars::VM)].setValue(0.86);
        regca.y()[index(Vars::IQ)].setValue(-0.2);
        regca.y()[index(Vars::IP)].setValue(0.85);
        regca.y()[index(Vars::VT)].setValue(0.98);
        regca.y()[index(Vars::II)].setValue(0.18);
        regca.y()[index(Vars::IQEXTRA)].setValue(0.03);
        regca.y()[index(Vars::IL)].setValue(0.72);
        regca.y()[index(Vars::IR)].setValue(0.5);
        regca.y()[index(Vars::LP)].setValue(-0.4);
        regca.y()[index(Vars::UP)].setValue(0.3);
        regca.y()[index(Vars::PBR)].setValue(0.52);
        regca.y()[index(Vars::QBR)].setValue(-0.046);

        regca.yp()[index(Vars::VM)].setValue(0.01);
        regca.yp()[index(Vars::IQ)].setValue(-0.02);
        regca.yp()[index(Vars::IP)].setValue(0.03);
      }

      std::vector<DependencyTracking::Variable::DependencyMap> dependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        auto data = makeDynamicData();

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

        for (IdxT i = 0; i < regca.size(); ++i)
        {
          regca.y()[static_cast<size_t>(i)].setVariableNumber(i);
          regca.yp()[static_cast<size_t>(i)].setVariableNumber(i);
        }
        for (IdxT i = 0; i < bus.size(); ++i)
        {
          bus.y()[static_cast<size_t>(i)].setVariableNumber(i + regca.size());
        }
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
        for (IdxT i = 0; i < regca.size(); ++i)
        {
          dependencies[static_cast<size_t>(i)] =
              regca.getResidual()[static_cast<size_t>(i)].getDependencies();
        }
        dependencies[static_cast<size_t>(regca.size())]     = bus.Ir().getDependencies();
        dependencies[static_cast<size_t>(regca.size() + 1)] = bus.Ii().getDependencies();

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> enzymeJacobian()
      {
        auto data = makeDynamicData();

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

        auto model_jacobian = regca.getJacobian();
        model_jacobian.deduplicate();
        return MapFromCOO(model_jacobian);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
