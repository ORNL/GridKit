#pragma once

#include <algorithm>
#include <iostream>
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
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ConverterRegcaTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterRegcaTests()  = default;
      ~ConverterRegcaTests() = default;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-10);

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> minimal(&bus);
        success *= (minimal.size() == static_cast<IdxT>(PhasorDynamics::Converter::RegcaInternalVariables::MAXIMUM));
        success *= (minimal.getMonitor() == nullptr);
        success *= (minimal.verify() > 0);

        auto                                            data = makeTestData();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> from_data(&bus, data);
        success *= (from_data.size() == static_cast<IdxT>(PhasorDynamics::Converter::RegcaInternalVariables::MAXIMUM));
        success *= (from_data.getMonitor() != nullptr);
        success *= (from_data.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        auto missing = makeTestData();
        missing.parameters.erase(Params::Tg);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> missing_model(&bus, missing);
        success *= (missing_model.verify() > 0);

        auto bad_switch                   = makeTestData();
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

      TestOutcome lifecycle()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        auto                                            data = makeTestData();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        success *= (regca.allocate() == 0);
        success *= (regca.initialize() == 0);
        success *= (regca.tagDifferentiable() == 0);
        success *= (regca.evaluateResidual() == 0);
        success *= (regca.evaluateJacobian() == 0);

        for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
        {
          success *= isEqual(regca.yp()[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(regca.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
        }

        using Vars = PhasorDynamics::Converter::RegcaInternalVariables;
        for (size_t i = 0; i < static_cast<size_t>(Vars::MAXIMUM); ++i)
        {
          const bool expected = (i == static_cast<size_t>(Vars::VM))
                                || (i == static_cast<size_t>(Vars::IQ))
                                || (i == static_cast<size_t>(Vars::IP));
          success *= (regca.tag()[i] == expected);
        }

        success *= isEqual(bus.Ir(), static_cast<ScalarT>(0.0), kTol);
        success *= isEqual(bus.Ii(), static_cast<ScalarT>(0.0), kTol);

        return success.report(__func__);
      }

      TestOutcome steadyStateInitializationGolden()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        auto data = makeGoldenTestData(true);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.92};
        ScalarT                                   iqcmd_value{-0.44};
        IdxT                                      ipcmd_index = 21;
        IdxT                                      iqcmd_index = 22;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        regca.allocate();
        success *= (regca.initialize() == 0);

        std::vector<ScalarT> expected_y(static_cast<size_t>(Vars::MAXIMUM), static_cast<ScalarT>(0.0));
        expected_y[index(Vars::VM)]      = static_cast<ScalarT>(1.0);
        expected_y[index(Vars::IQ)]      = static_cast<ScalarT>(-0.44);
        expected_y[index(Vars::IP)]      = static_cast<ScalarT>(0.92);
        expected_y[index(Vars::VT)]      = static_cast<ScalarT>(1.0);
        expected_y[index(Vars::II)]      = static_cast<ScalarT>(0.904);
        expected_y[index(Vars::IQEXTRA)] = static_cast<ScalarT>(0.0);
        expected_y[index(Vars::IL)]      = static_cast<ScalarT>(1.1);
        expected_y[index(Vars::IR)]      = static_cast<ScalarT>(0.472);
        expected_y[index(Vars::LP)]      = static_cast<ScalarT>(-70.0);
        expected_y[index(Vars::UP)]      = static_cast<ScalarT>(0.7);
        expected_y[index(Vars::PBR)]     = static_cast<ScalarT>(0.92);
        expected_y[index(Vars::QBR)]     = static_cast<ScalarT>(-0.44);

        success *= vectorMatches(regca.y(), expected_y, "REGCA initialization state");
        for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
        {
          success *= isEqual(regca.yp()[i], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }

      TestOutcome attachedSignalInitialization()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);
        bus.allocate();
        bus.initialize();

        auto                                            data = makeTestData();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.75};
        ScalarT                                   iqcmd_value{-0.20};
        IdxT                                      ipcmd_index = 21;
        IdxT                                      iqcmd_index = 22;

        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        regca.allocate();
        success *= (regca.initialize() == 0);

        success *= isEqual(ipcmd_value, static_cast<ScalarT>(0.75), kTol);
        success *= isEqual(iqcmd_value, static_cast<ScalarT>(-0.20), kTol);
        success *= isEqual(regca.y()[index(Vars::IP)], static_cast<ScalarT>(0.75), kTol);
        success *= isEqual(regca.y()[index(Vars::IQ)], static_cast<ScalarT>(-0.20), kTol);

        return success.report(__func__);
      }

      TestOutcome fallbackInit()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        auto data                    = makeTestData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);
        data.parameters[Params::IL1] = static_cast<RealT>(2.0);
        data.parameters[Params::P0]  = static_cast<RealT>(0.8);
        data.parameters[Params::Q0]  = static_cast<RealT>(-0.1);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);
        regca.allocate();
        success *= (regca.initialize() == 0);

        success *= isEqual(regca.y()[index(Vars::IR)], static_cast<ScalarT>(1.16), kTol);
        success *= isEqual(regca.y()[index(Vars::II)], static_cast<ScalarT>(1.12), kTol);
        success *= isEqual(regca.y()[index(Vars::PBR)], static_cast<ScalarT>(0.8), kTol);
        success *= isEqual(regca.y()[index(Vars::QBR)], static_cast<ScalarT>(-0.1), kTol);
        success *= isEqual(regca.y()[index(Vars::IP)], static_cast<ScalarT>(1.6), kTol);
        success *= isEqual(regca.y()[index(Vars::IQ)], static_cast<ScalarT>(-0.2), kTol);

        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);

        for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
        {
          success *= isEqual(regca.getResidual()[i], static_cast<ScalarT>(0.0), kTol);
        }
        success *= isEqual(bus.Ir(), static_cast<ScalarT>(0.58), kTol);
        success *= isEqual(bus.Ii(), static_cast<ScalarT>(0.56), kTol);

        return success.report(__func__);
      }

      TestOutcome invalidInitialization()
      {
        TestStatus success = true;

        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT> high_voltage_bus(1.3, 0.0);
        high_voltage_bus.allocate();
        high_voltage_bus.initialize();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> high_voltage_regca(&high_voltage_bus, data);
        high_voltage_regca.allocate();
        success *= (high_voltage_regca.initialize() > 0);

        auto low_voltage_fallback_data                   = data;
        low_voltage_fallback_data.parameters[Params::P0] = static_cast<RealT>(0.75);
        PhasorDynamics::Bus<ScalarT, IdxT> low_voltage_fallback_bus(0.2, 0.0);
        low_voltage_fallback_bus.allocate();
        low_voltage_fallback_bus.initialize();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> low_voltage_fallback_regca(
            &low_voltage_fallback_bus, low_voltage_fallback_data);
        low_voltage_fallback_regca.allocate();
        success *= (low_voltage_fallback_regca.initialize() > 0);

        PhasorDynamics::Bus<ScalarT, IdxT> low_voltage_attached_bus(0.2, 0.0);
        low_voltage_attached_bus.allocate();
        low_voltage_attached_bus.initialize();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> low_voltage_attached_regca(
            &low_voltage_attached_bus, low_voltage_fallback_data);
        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        ScalarT                                   ipcmd_value{0.75};
        IdxT                                      ipcmd_index = 21;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        low_voltage_attached_regca.getSignals().template attachSignalNode<Ext::IPCMD>(
            &ipcmd_node);
        low_voltage_attached_regca.allocate();
        success *= (low_voltage_attached_regca.initialize() == 0);

        return success.report(__func__);
      }

      TestOutcome signalVerification()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, makeTestData());

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.25};
        ScalarT                                   iqcmd_value{-0.10};
        IdxT                                      ipcmd_index = 21;
        IdxT                                      iqcmd_index = 22;

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

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(nullptr, makeTestData());
        success *= (regca.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome residualGoldenVectors()
      {
        TestStatus success = true;

        std::vector<ScalarT> positive_q_lvpl(static_cast<size_t>(Vars::MAXIMUM), static_cast<ScalarT>(0.0));
        positive_q_lvpl[index(Vars::VM)]      = static_cast<ScalarT>(0.29);
        positive_q_lvpl[index(Vars::IQ)]      = static_cast<ScalarT>(0.52);
        positive_q_lvpl[index(Vars::IP)]      = static_cast<ScalarT>(0.21999997439919405);
        positive_q_lvpl[index(Vars::VT)]      = static_cast<ScalarT>(0.0046000000000000485);
        positive_q_lvpl[index(Vars::II)]      = static_cast<ScalarT>(0.2546);
        positive_q_lvpl[index(Vars::IQEXTRA)] = static_cast<ScalarT>(-0.03);
        positive_q_lvpl[index(Vars::IL)]      = static_cast<ScalarT>(0.29199937917427254);
        positive_q_lvpl[index(Vars::IR)]      = static_cast<ScalarT>(0.26);
        positive_q_lvpl[index(Vars::LP)]      = static_cast<ScalarT>(-69.6);
        positive_q_lvpl[index(Vars::UP)]      = static_cast<ScalarT>(-0.2999999999999803);
        positive_q_lvpl[index(Vars::PBR)]     = static_cast<ScalarT>(0.0);
        positive_q_lvpl[index(Vars::QBR)]     = static_cast<ScalarT>(0.0);

        std::vector<ScalarT> nonpositive_q_no_lvpl(static_cast<size_t>(Vars::MAXIMUM), static_cast<ScalarT>(0.0));
        nonpositive_q_no_lvpl[index(Vars::VM)]      = static_cast<ScalarT>(0.29);
        nonpositive_q_no_lvpl[index(Vars::IQ)]      = static_cast<ScalarT>(1.5200000000000002);
        nonpositive_q_no_lvpl[index(Vars::IP)]      = static_cast<ScalarT>(0.21999997439919405);
        nonpositive_q_no_lvpl[index(Vars::VT)]      = static_cast<ScalarT>(0.0046000000000000485);
        nonpositive_q_no_lvpl[index(Vars::II)]      = static_cast<ScalarT>(0.2546);
        nonpositive_q_no_lvpl[index(Vars::IQEXTRA)] = static_cast<ScalarT>(-0.03);
        nonpositive_q_no_lvpl[index(Vars::IL)]      = static_cast<ScalarT>(0.29199937917427254);
        nonpositive_q_no_lvpl[index(Vars::IR)]      = static_cast<ScalarT>(0.26);
        nonpositive_q_no_lvpl[index(Vars::LP)]      = static_cast<ScalarT>(-69.6);
        nonpositive_q_no_lvpl[index(Vars::UP)]      = static_cast<ScalarT>(0.39999999999999997);
        nonpositive_q_no_lvpl[index(Vars::PBR)]     = static_cast<ScalarT>(0.0);
        nonpositive_q_no_lvpl[index(Vars::QBR)]     = static_cast<ScalarT>(0.0);

        success *= residualGoldenVectorCase(static_cast<ScalarT>(0.1), true, positive_q_lvpl);
        success *= residualGoldenVectorCase(static_cast<ScalarT>(-0.1), false, nonpositive_q_no_lvpl);

        return success.report(__func__);
      }

      TestOutcome busInjection()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.8, 0.6);
        bus.allocate();
        bus.initialize();

        auto data                    = makeTestData();
        data.parameters[Params::mva] = static_cast<RealT>(50.0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> pbranch_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> qbranch_node;
        ScalarT                                   ipcmd_value{1.6};
        ScalarT                                   iqcmd_value{-0.2};
        IdxT                                      ipcmd_index = 21;
        IdxT                                      iqcmd_index = 22;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);
        regca.getSignals().template assignSignalNode<Vars::PBR>(&pbranch_node);
        regca.getSignals().template assignSignalNode<Vars::QBR>(&qbranch_node);

        regca.allocate();
        success *= (regca.initialize() == 0);

        bus.evaluateResidual();
        success *= (regca.evaluateResidual() == 0);

        success *= isEqual(bus.Ir(), static_cast<ScalarT>(0.58), kTol);
        success *= isEqual(bus.Ii(), static_cast<ScalarT>(0.56), kTol);
        success *= isEqual(pbranch_node.read(), static_cast<ScalarT>(0.8), kTol);
        success *= isEqual(qbranch_node.read(), static_cast<ScalarT>(-0.1), kTol);

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
    { "number": 1, "class": "bus", "name": "Bus 1", "init": { "Vr": 1.0, "Vi": 0.0 }, "v_base": 1.0 }
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
        success                *= (regca_data.ports.at(PhasorDynamics::Converter::RegcaPorts::bus) == 1);
        success                *= (std::get_if<double>(&regca_data.parameters.at(Params::P0)) != nullptr);
        success                *= (std::get_if<double>(&regca_data.parameters.at(Params::Q0)) != nullptr);
        success                *= (std::get_if<size_t>(&regca_data.parameters.at(Params::mva)) != nullptr);
        success                *= (std::get_if<bool>(&regca_data.parameters.at(Params::sL)) != nullptr);

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

        auto dependency_tracking_jacobian = DependencyTrackingJacobian();
        auto enzyme_jacobian              = EnzymeJacobian();

        success          *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
        const auto nrows  = std::min(dependency_tracking_jacobian.size(), enzyme_jacobian.size());
        for (size_t i = 0; i < nrows; ++i)
        {
          success *= isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], static_cast<RealT>(1.0e-8));
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

      auto makeTestData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        using Ports = PhasorDynamics::Converter::RegcaPorts;
        using Mon   = PhasorDynamics::Converter::RegcaMonitorableVariables;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> data;
        data.device_class          = "Regca";
        data.disambiguation_string = "regca_test";
        data.ports[Ports::bus]     = 1;
        data.monitored_variables.insert(Mon::ir);
        data.monitored_variables.insert(Mon::ii);

        data.parameters[Params::mva]    = static_cast<IdxT>(100);
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

      auto makeGoldenTestData(bool use_lvpl) -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        auto data                       = makeTestData();
        data.parameters[Params::Tg]     = static_cast<RealT>(0.2);
        data.parameters[Params::TM]     = static_cast<RealT>(0.4);
        data.parameters[Params::Rqmax]  = static_cast<RealT>(0.5);
        data.parameters[Params::Rqmin]  = static_cast<RealT>(-0.6);
        data.parameters[Params::Rpmax]  = static_cast<RealT>(0.7);
        data.parameters[Params::sL]     = use_lvpl;
        data.parameters[Params::IL1]    = static_cast<RealT>(1.1);
        data.parameters[Params::VL0]    = static_cast<RealT>(0.4);
        data.parameters[Params::VL1]    = static_cast<RealT>(0.9);
        data.parameters[Params::VA0]    = static_cast<RealT>(0.4);
        data.parameters[Params::VA1]    = static_cast<RealT>(0.9);
        data.parameters[Params::Vhvmax] = static_cast<RealT>(1.3);
        return data;
      }

      bool invalidParameterCase(PhasorDynamics::Bus<ScalarT, IdxT>& bus, Params param, RealT value)
      {
        auto data              = makeTestData();
        data.parameters[param] = value;
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> model(&bus, data);
        return model.verify() > 0;
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

      void setGoldenResidualState(PhasorDynamics::Converter::Regca<ScalarT, IdxT>& regca)
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

      bool residualGoldenVectorCase(ScalarT initial_iqcmd, bool use_lvpl, const std::vector<ScalarT>& expected)
      {
        bool success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(0.95, 0.25);
        bus.allocate();
        bus.initialize();

        auto                                            data = makeGoldenTestData(use_lvpl);
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(&bus, data);

        PhasorDynamics::SignalNode<ScalarT, IdxT> ipcmd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> iqcmd_node;
        ScalarT                                   ipcmd_value{0.9};
        ScalarT                                   iqcmd_value{initial_iqcmd};
        IdxT                                      ipcmd_index = 21;
        IdxT                                      iqcmd_index = 22;
        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        regca.getSignals().template attachSignalNode<Ext::IPCMD>(&ipcmd_node);
        regca.getSignals().template attachSignalNode<Ext::IQCMD>(&iqcmd_node);

        regca.allocate();
        regca.initialize();
        iqcmd_value = static_cast<ScalarT>(0.1);
        setGoldenResidualState(regca);

        bus.evaluateResidual();
        regca.evaluateResidual();

        success *= vectorMatches(regca.getResidual(), expected, "REGCA residual");

        return success;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      void setJacobianState(PhasorDynamics::Converter::Regca<ScalarT, IdxT>& regca,
                            PhasorDynamics::Bus<ScalarT, IdxT>&              bus)
      {
        bus.y()[0] = static_cast<ScalarT>(0.95);
        bus.y()[1] = static_cast<ScalarT>(0.25);

        setGoldenResidualState(regca);
      }

      void setJacobianStateDep(
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

      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        using DepVar = DependencyTracking::Variable;

        auto data = makeGoldenTestData(true);

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
        setJacobianStateDep(regca, bus);

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
          dependencies[static_cast<size_t>(i)] = regca.getResidual()[static_cast<size_t>(i)].getDependencies();
        }
        dependencies[static_cast<size_t>(regca.size())]     = bus.Ir().getDependencies();
        dependencies[static_cast<size_t>(regca.size() + 1)] = bus.Ii().getDependencies();

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        auto data = makeGoldenTestData(true);

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
        setJacobianState(regca, bus);
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
