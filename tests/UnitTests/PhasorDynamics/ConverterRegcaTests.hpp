#pragma once

#include <iostream>
#include <limits>
#include <sstream>
#include <variant>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
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
    class ConverterRegcaTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ConverterRegcaTests()  = default;
      ~ConverterRegcaTests() = default;

      static constexpr ScalarT kTol = static_cast<ScalarT>(1.0e-14);

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> minimal(&bus);
        success *= (minimal.size() == static_cast<IdxT>(PhasorDynamics::Converter::RegcaInternalVariables::MAXIMUM));
        success *= (minimal.getMonitor() == nullptr);

        auto                                            data = makeTestData();
        PhasorDynamics::Converter::Regca<ScalarT, IdxT> from_data(&bus, data);
        success *= (from_data.size() == static_cast<IdxT>(PhasorDynamics::Converter::RegcaInternalVariables::MAXIMUM));
        success *= (from_data.getMonitor() != nullptr);

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

        const auto& y  = regca.y();
        const auto& yp = regca.yp();
        const auto& f  = regca.getResidual();
        for (size_t i = 0; i < static_cast<size_t>(regca.size()); ++i)
        {
          success *= isEqual(y[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(yp[i], static_cast<ScalarT>(0.0), kTol);
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), kTol);
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

        regca.getSignals().template attachSignalNode<PhasorDynamics::Converter::RegcaExternalVariables::IPCMD>(&ipcmd_node);
        success *= (regca.verify() > 0);

        ipcmd_node.set(&ipcmd_value, &ipcmd_index);
        success *= (regca.verify() == 0);

        regca.getSignals().template attachSignalNode<PhasorDynamics::Converter::RegcaExternalVariables::IQCMD>(&iqcmd_node);
        success *= (regca.verify() > 0);

        iqcmd_node.set(&iqcmd_value, &iqcmd_index);
        success *= (regca.verify() == 0);

        return success.report(__func__);
      }

      TestOutcome nullBusVerification()
      {
        TestStatus success = true;

        PhasorDynamics::Converter::Regca<ScalarT, IdxT> regca(nullptr);
        success *= (regca.verify() > 0);

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
    "case_name": "REGCA skeleton",
    "case_description": "REGCA parser smoke test",
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
        "P0": 1.0,
        "Q0": 0.0,
        "Sconv": 100,
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
      "mon": ["ir", "ii"]
    }
  ]
}
)json");

        auto data  = PhasorDynamics::parseSystemModelData(input);
        success   *= (data.regca.size() == 1);
        success   *= (data.regca[0].device_class == "Regca");
        success   *= (data.regca[0].ports.at(PhasorDynamics::Converter::RegcaPorts::bus) == 1);
        success   *= (std::get_if<size_t>(&data.regca[0].parameters.at(PhasorDynamics::Converter::RegcaParameters::Sconv)) != nullptr);
        success   *= (std::get_if<bool>(&data.regca[0].parameters.at(PhasorDynamics::Converter::RegcaParameters::sL)) != nullptr);

        PhasorDynamics::SystemModel<double, size_t> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.tagDifferentiable() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.evaluateJacobian() == 0);
        success *= (system.size() == 12);
        success *= isEqual(system.getResidual()[0], 0.0, static_cast<double>(kTol));
        success *= isEqual(system.getResidual()[1], 0.0, static_cast<double>(kTol));

        return success.report(__func__);
      }

    private:
      auto makeTestData() -> PhasorDynamics::Converter::RegcaData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Converter::RegcaParameters;
        using Ports  = PhasorDynamics::Converter::RegcaPorts;
        using Mon    = PhasorDynamics::Converter::RegcaMonitorableVariables;

        PhasorDynamics::Converter::RegcaData<RealT, IdxT> data;
        data.device_class          = "Regca";
        data.disambiguation_string = "regca_test";
        data.ports[Ports::bus]     = 1;
        data.monitored_variables.insert(Mon::ir);
        data.monitored_variables.insert(Mon::ii);

        data.parameters[Params::P0]     = static_cast<RealT>(1.0);
        data.parameters[Params::Q0]     = static_cast<RealT>(0.0);
        data.parameters[Params::Sconv]  = static_cast<IdxT>(100);
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
    };
  } // namespace Testing
} // namespace GridKit
