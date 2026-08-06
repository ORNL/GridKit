#pragma once

#include <iostream>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/GensalData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    using json = nlohmann::json;

    template <typename RealT, typename IdxT>
    class CaseFormatTests
    {
    public:
      CaseFormatTests()  = default;
      ~CaseFormatTests() = default;

      TestOutcome simpleParse()
      {
        using namespace GridKit::PhasorDynamics;
        using namespace GridKit::Model;
        using BusData   = BusData<RealT, IdxT>;
        using BusType   = typename BusData::BusType;
        using RegcaData = Converter::RegcaData<RealT, IdxT>;

        const char data[] =
            R"({
               "header": {
                   "format_version": 0.2,
                   "format_revision": 0,
                   "case_name": "Two-bus test case 1",
                   "case_description": "A two-bus test case for demonstrating the dynamics format",
                   "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed"
               },
               "params": {
                   "freq_base": 50.0,
                   "va_base": 50000000.0
               },
               "monitors": [
                   {
                       "file_name": "mon.csv",
                       "format": "csv"
                   },
                   {
                       "file_name": "mon.json",
                       "format": "json"
                   }
               ],
               "buses": [
                   { "number": 1, "class": "Bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "params": {"kv": 115.0}, "mon": ["Vr", "Vi"] },
                   { "number": 2, "class": "BusInfinite", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0.0}, "params": {"kv": 115.0} }
               ],
               "devices": [
                   { "class": "Branch", "ports": {"bus1":1, "bus2":2}, "id": "1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0, "Gmag":0.01, "Bmag":0.02, "tap":1.05, "phase":0.1} },
                   { "class": "Genrou", "ports": {"bus":1}, "id": "1", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05,
                          "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                   { "class": "Gensal", "ports": {"bus":1}, "id": "2", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05,
                          "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                   { "class": "BusFault", "ports": {"bus":1}, "id": "1", "params": {"state0": false, "R":0.0, "X":1e-3} },
                   { "class": "Regca", "ports": {"bus":1}, "id": "CV1", "params": {"p0":0.0, "q0":0.0, "mva":100, "Tg":0.02, "TM":0.02, "Rqmax":999.0, "Rqmin":-999.0, "Rpmax":999.0, "sL":true, "IL1":1.1, "VL0":0.4, "VL1":0.9, "VA0":0.4, "VA1":0.9, "Vhvmax":1.2}, "mon": ["ir", "ii", "p", "q"] }
               ]
            })";

        TestStatus                   success = true;
        SystemModelData<RealT, IdxT> result  = json::parse(data);

        success *= result.format_version == 0.2;
        success *= result.format_revision == 0;
        success *= result.case_name == "Two-bus test case 1";
        success *= result.case_date_time == std::nullopt;
        success *= result.case_description == "A two-bus test case for demonstrating the dynamics format";
        success *= result.case_comments == "This case is set up to monitor the voltage at both buses and the machine angle and speed";
        success *= result.freq_base == 50.0;
        success *= result.va_base == 50.0e6;

        success *= result.monitor_sink[0].file_name == "mon.csv";
        success *= result.monitor_sink[0].format == VariableMonitorFormat::CSV;
        success *= result.monitor_sink[1].file_name == "mon.json";
        success *= result.monitor_sink[1].format == VariableMonitorFormat::JSON;

        success *= result.bus.size() == 2;
        success *= result.branch.size() == 1;
        success *= result.bus_fault.size() == 1;
        success *= result.genrou.size() == 1;
        success *= result.gensal.size() == 1;
        success *= result.regca.size() == 1;
        success *= result.loadz.size() == 0;

        success *= result.bus[0].bus_id == 1;
        success *= result.bus[0].bus_type == BusType::DEFAULT;
        success *= result.bus[0].name == "Bus 1";
        success *= result.bus[0].Vr0 == 0.994988;
        success *= result.bus[0].Vi0 == 0.099997;
        success *= std::get<RealT>(result.bus[0].parameters[BusData::Parameters::kv]) == 115.0;
        success *= result.bus[0].monitored_variables.contains(BusData::MonitorableVariables::Vr);
        success *= result.bus[0].monitored_variables.contains(BusData::MonitorableVariables::Vi);
        success *= result.bus[1].bus_id == 2;
        success *= result.bus[1].bus_type == BusType::SLACK;
        success *= result.bus[1].name == "Bus 2";
        success *= result.bus[1].Vr0 == 1.0;
        success *= result.bus[1].Vi0 == 0.0;
        success *= std::get<RealT>(result.bus[1].parameters[BusData::Parameters::kv]) == 115.0;
        success *= result.bus[1].monitored_variables.empty();

        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::R]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::X]) == 0.1;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::G]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::B]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::Gmag]) == 0.01;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::Bmag]) == 0.02;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::tap]) == 1.05;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::phase]) == 0.1;
        success *= result.branch[0].buses[BranchBuses::bus1] == 1;
        success *= result.branch[0].buses[BranchBuses::bus2] == 2;
        success *= result.branch[0].disambiguation_string == "1";
        success *= result.branch[0].monitored_variables.empty();

        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::p0]) == 1.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::q0]) == 0.05013;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::H]) == 3.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::D]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Ra]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tdop]) == 7.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tdopp]) == 0.04;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tqop]) == 0.75;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tqopp]) == 0.05;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xd]) == 2.1;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xdp]) == 0.2;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xdpp]) == 0.18;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xq]) == 0.5;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xqp]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xqpp]) == 0.18;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xl]) == 0.15;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::S10]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::S12]) == 0.0;
        success *= result.genrou[0].buses[GenrouBuses::bus] == 1;
        success *= result.genrou[0].disambiguation_string == "1";
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::delta);
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::omega);

        success *= std::get<RealT>(result.gensal[0].parameters[GensalParameters::p0]) == 1.0;
        success *= std::get<RealT>(result.gensal[0].parameters[GensalParameters::q0]) == 0.05013;
        success *= std::get<RealT>(result.gensal[0].parameters[GensalParameters::Xd]) == 2.1;
        success *= std::get<RealT>(result.gensal[0].parameters[GensalParameters::Xq]) == 0.5;
        success *= result.gensal[0].buses[GensalBuses::bus] == 1;
        success *= result.gensal[0].disambiguation_string == "2";
        success *= result.gensal[0].monitored_variables.contains(GensalMonitorableVariables::delta);
        success *= result.gensal[0].monitored_variables.contains(GensalMonitorableVariables::omega);

        success *= std::get<RealT>(result.regca[0].parameters[RegcaData::Parameters::p0]) == 0.0;
        success *= std::get<RealT>(result.regca[0].parameters[RegcaData::Parameters::q0]) == 0.0;
        success *= std::get<IdxT>(result.regca[0].parameters[RegcaData::Parameters::mva]) == 100;
        success *= std::get<bool>(result.regca[0].parameters[RegcaData::Parameters::sL]);
        success *= result.regca[0].buses[RegcaData::Buses::bus] == 1;
        success *= result.regca[0].disambiguation_string == "CV1";
        success *= result.regca[0].signal_inputs.empty();
        success *= result.regca[0].signal_outputs.empty();
        success *= result.regca[0].monitored_variables.contains(RegcaData::MonitorableVariables::ir);
        success *= result.regca[0].monitored_variables.contains(RegcaData::MonitorableVariables::ii);
        success *= result.regca[0].monitored_variables.contains(RegcaData::MonitorableVariables::p);
        success *= result.regca[0].monitored_variables.contains(RegcaData::MonitorableVariables::q);

        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::R]) == 0.0;
        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::X]) == 1e-3;
        success *= !std::get<bool>(result.bus_fault[0].parameters[BusFaultParameters::state0]);
        success *= result.bus_fault[0].buses[BusFaultBuses::bus] == 1;
        success *= result.bus_fault[0].disambiguation_string == "1";
        success *= result.bus_fault[0].monitored_variables.empty();

        return success.report(__func__);
      }

      TestOutcome signalParse()
      {
        using namespace GridKit::PhasorDynamics;
        using BusData    = BusData<RealT, IdxT>;
        using BusType    = typename BusData::BusType;
        using Esdc1aData = Exciter::Esdc1aData<RealT, IdxT>;

        const char data[] =
            R"({
               "header": {
                   "format_version": 0.2,
                   "format_revision": 0,
                   "case_name": "Two-bus test case 2",
                   "case_description": "A two-bus test case for demonstrating the dynamics format",
                   "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed"
               },
               "params": {
                   "freq_base": 60.0,
                   "va_base": 100000000.0
               },
               "buses": [
                   { "number": 1, "class": "Bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "params": {"kv": 115.0}, "mon": ["Vr", "Vi"] },
                   { "number": 2, "class": "BusInfinite", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0.0}, "params": {"kv": 115.0} }
               ],
               "signals": [
                   { "signal_id": 1, "name": "Machine Speed Deviation"},
                   { "signal_id": 2, "name": "Mechanical Power"},
                   { "signal_id": 3, "name": "Excitation Field"},
                   { "signal_id": 4, "name": "Voltage Reference"},
                   { "signal_id": 5, "name": "Stabilizer Signal"},
                   { "signal_id": 6, "name": "Under-excitation Limiter"}
               ],
               "devices": [
                   { "class": "Branch", "ports": {"bus1":1, "bus2":2}, "id": "BR1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0, "tap":1.05, "phase":0.1} },
                   { "class": "Genrou", "ports": {"bus":1, "speed": 1, "pmech":2, "efd":3}, "id": "DV1", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05, "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                   { "class": "Tgov1", "ports": {"speed": 1, "pmech":2}, "id": "DV2", "params": {"R":0.05, "T1":0.5,"T2":2.5, "T3":7.5, "Pvmax":0.0, "Pvmin":1.0, "Dt":0.0}},
                   { "class": "Esdc1a", "ports": {"bus":1, "speed":1, "vref":4, "vs":5, "vuel":6, "efd":3}, "id": "DV5", "params": {"Tr":0.0, "Ka":40.0, "Ta":0.1, "Tb":0.0, "Tc":0.0, "Vrmax":1.0, "Vrmin":-1.0, "Ke":0.1, "Te":0.5, "Kf":0.05, "Tf1":0.7, "Spdmlt":false, "E1":2.8, "Se1":0.08, "E2":3.7, "Se2":0.33, "UEL":0, "exclim":true}, "mon": ["efd", "vc", "vr", "vf", "se", "vfe"] },
                   { "class": "Ieeet1", "ports": {"bus":1, "speed": 1, "efd":3}, "id": "DV3", "params": {"Tr":0.0, "Ka":50.0, "Ta":0.04, "Ke":-0.06, "Te":0.6, "Kf":0.09, "Tf":1.46, "Vrmin":-1.0, "Vrmax":1.0, "E1":2.8, "E2":3.373, "Se1":0.04, "Se2":0.33, "Ispdlim":0.0}},
                   { "class": "SexsPti", "ports": {"bus":1, "efd":3}, "id": "DV4", "params": {"Ta":0.1, "Tb":0.5, "Te":0.8, "K":10.0, "Efdmax":5.0, "Efdmin":-5.0}},
                   { "class": "BusFault", "ports": {"bus":1}, "id": "1", "params": {"state0": false, "R":0.0, "X":1e-3} }
               ]
            })";

        TestStatus                   success = true;
        SystemModelData<RealT, IdxT> result  = json::parse(data);

        success *= result.format_version == 0.2;
        success *= result.format_revision == 0;
        success *= result.case_name == "Two-bus test case 2";
        success *= result.case_date_time == std::nullopt;
        success *= result.case_description == "A two-bus test case for demonstrating the dynamics format";
        success *= result.case_comments == "This case is set up to monitor the voltage at both buses and the machine angle and speed";
        success *= result.freq_base == 60.0;
        success *= result.va_base == 100.0e6;

        success *= result.bus.size() == 2;
        success *= result.branch.size() == 1;
        success *= result.bus_fault.size() == 1;
        success *= result.genrou.size() == 1;
        success *= result.gov.size() == 1;
        success *= result.esdc1a.size() == 1;
        success *= result.loadz.size() == 0;
        success *= result.exciter.size() == 1;
        success *= result.sexspti.size() == 1;

        success *= result.bus[0].bus_id == 1;
        success *= result.bus[0].bus_type == BusType::DEFAULT;
        success *= result.bus[0].name == "Bus 1";
        success *= result.bus[0].Vr0 == 0.994988;
        success *= result.bus[0].Vi0 == 0.099997;
        success *= std::get<RealT>(result.bus[0].parameters[BusData::Parameters::kv]) == 115.0;
        success *= result.bus[0].monitored_variables.contains(BusData::MonitorableVariables::Vr);
        success *= result.bus[0].monitored_variables.contains(BusData::MonitorableVariables::Vi);
        success *= result.bus[1].bus_id == 2;
        success *= result.bus[1].bus_type == BusType::SLACK;
        success *= result.bus[1].name == "Bus 2";
        success *= result.bus[1].Vr0 == 1.0;
        success *= result.bus[1].Vi0 == 0.0;
        success *= std::get<RealT>(result.bus[1].parameters[BusData::Parameters::kv]) == 115.0;
        success *= result.bus[1].monitored_variables.empty();

        success *= result.signal[0].signal_id == 1;
        success *= result.signal[0].name == "Machine Speed Deviation";
        success *= result.signal[1].signal_id == 2;
        success *= result.signal[1].name == "Mechanical Power";
        success *= result.signal[2].signal_id == 3;
        success *= result.signal[2].name == "Excitation Field";
        success *= result.signal[3].signal_id == 4;
        success *= result.signal[3].name == "Voltage Reference";
        success *= result.signal[4].signal_id == 5;
        success *= result.signal[4].name == "Stabilizer Signal";
        success *= result.signal[5].signal_id == 6;
        success *= result.signal[5].name == "Under-excitation Limiter";

        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::R]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::X]) == 0.1;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::G]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::B]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::tap]) == 1.05;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::phase]) == 0.1;
        success *= result.branch[0].buses[BranchBuses::bus1] == 1;
        success *= result.branch[0].buses[BranchBuses::bus2] == 2;
        success *= result.branch[0].disambiguation_string == "BR1";
        success *= result.branch[0].monitored_variables.empty();

        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::p0]) == 1.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::q0]) == 0.05013;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::H]) == 3.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::D]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Ra]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tdop]) == 7.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tdopp]) == 0.04;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tqop]) == 0.75;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Tqopp]) == 0.05;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xd]) == 2.1;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xdp]) == 0.2;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xdpp]) == 0.18;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xq]) == 0.5;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xqp]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xqpp]) == 0.18;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::Xl]) == 0.15;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::S10]) == 0.0;
        success *= std::get<RealT>(result.genrou[0].parameters[GenrouParameters::S12]) == 0.0;
        success *= result.genrou[0].buses[GenrouBuses::bus] == 1;
        success *= result.genrou[0].signal_outputs[GenrouSignalOutputs::speed] == 1;
        success *= result.genrou[0].signal_inputs[GenrouSignalInputs::pmech] == 2;
        success *= result.genrou[0].signal_inputs[GenrouSignalInputs::efd] == 3;
        success *= result.genrou[0].disambiguation_string == "DV1";
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::delta);
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::omega);

        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::R]) == 0.05;
        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::T1]) == 0.5;
        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::T2]) == 2.5;
        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::T3]) == 7.5;
        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::Pvmax]) == 0;
        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::Pvmin]) == 1;
        success *= std::get<RealT>(result.gov[0].parameters[Governor::Tgov1Parameters::Dt]) == 0;
        success *= result.gov[0].signal_inputs[Governor::Tgov1SignalInputs::speed] == 1;
        success *= result.gov[0].signal_outputs[Governor::Tgov1SignalOutputs::pmech] == 2;
        success *= result.gov[0].disambiguation_string == "DV2";

        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Tr]) == 0.0;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Ka]) == 40.0;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Ta]) == 0.1;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Tb]) == 0.0;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Tc]) == 0.0;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Vrmax]) == 1.0;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Vrmin]) == -1.0;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Ke]) == 0.1;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Te]) == 0.5;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Kf]) == 0.05;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Tf1]) == 0.7;
        success *= !std::get<bool>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Spdmlt]);
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::E1]) == 2.8;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Se1]) == 0.08;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::E2]) == 3.7;
        success *= std::get<RealT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::Se2]) == 0.33;
        success *= std::get<IdxT>(result.esdc1a[0].parameters[Esdc1aData::Parameters::UEL]) == 0;
        success *= std::get<bool>(result.esdc1a[0].parameters[Esdc1aData::Parameters::exclim]);
        success *= result.esdc1a[0].buses[Esdc1aData::Buses::bus] == 1;
        success *= result.esdc1a[0].signal_inputs[Esdc1aData::SignalInputs::speed] == 1;
        success *= result.esdc1a[0].signal_inputs[Esdc1aData::SignalInputs::vref] == 4;
        success *= result.esdc1a[0].signal_inputs[Esdc1aData::SignalInputs::vs] == 5;
        success *= result.esdc1a[0].signal_inputs[Esdc1aData::SignalInputs::vuel] == 6;
        success *= result.esdc1a[0].signal_outputs[Esdc1aData::SignalOutputs::efd] == 3;
        success *= result.esdc1a[0].disambiguation_string == "DV5";
        success *= result.esdc1a[0].monitored_variables.contains(Esdc1aData::MonitorableVariables::efd);
        success *= result.esdc1a[0].monitored_variables.contains(Esdc1aData::MonitorableVariables::vc);
        success *= result.esdc1a[0].monitored_variables.contains(Esdc1aData::MonitorableVariables::vr);
        success *= result.esdc1a[0].monitored_variables.contains(Esdc1aData::MonitorableVariables::vf);
        success *= result.esdc1a[0].monitored_variables.contains(Esdc1aData::MonitorableVariables::se);
        success *= result.esdc1a[0].monitored_variables.contains(Esdc1aData::MonitorableVariables::vfe);

        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Tr]) == 0.0;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Ka]) == 50.0;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Ta]) == 0.04;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Ke]) == -0.06;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Te]) == 0.6;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Kf]) == 0.09;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Tf]) == 1.46;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Vrmin]) == -1.0;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Vrmax]) == 1.0;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::E1]) == 2.8;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::E2]) == 3.373;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Se1]) == 0.04;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Se2]) == 0.33;
        success *= std::get<RealT>(result.exciter[0].parameters[Exciter::Ieeet1Parameters::Ispdlim]) == 0.0;
        success *= result.exciter[0].buses[Exciter::Ieeet1Buses::bus] == 1;
        success *= result.exciter[0].signal_inputs[Exciter::Ieeet1SignalInputs::speed] == 1;
        success *= result.exciter[0].signal_outputs[Exciter::Ieeet1SignalOutputs::efd] == 3;
        success *= result.exciter[0].disambiguation_string == "DV3";

        success *= std::get<RealT>(result.sexspti[0].parameters[Exciter::SexsPtiParameters::Ta]) == 0.1;
        success *= std::get<RealT>(result.sexspti[0].parameters[Exciter::SexsPtiParameters::Tb]) == 0.5;
        success *= std::get<RealT>(result.sexspti[0].parameters[Exciter::SexsPtiParameters::Te]) == 0.8;
        success *= std::get<RealT>(result.sexspti[0].parameters[Exciter::SexsPtiParameters::K]) == 10.0;
        success *= std::get<RealT>(result.sexspti[0].parameters[Exciter::SexsPtiParameters::Efdmax]) == 5.0;
        success *= std::get<RealT>(result.sexspti[0].parameters[Exciter::SexsPtiParameters::Efdmin]) == -5.0;
        success *= result.sexspti[0].buses[Exciter::SexsPtiBuses::bus] == 1;
        success *= result.sexspti[0].signal_outputs[Exciter::SexsPtiSignalOutputs::efd] == 3;
        success *= result.sexspti[0].disambiguation_string == "DV4";

        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::R]) == 0.0;
        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::X]) == 1e-3;
        success *= !std::get<bool>(result.bus_fault[0].parameters[BusFaultParameters::state0]);
        success *= result.bus_fault[0].buses[BusFaultBuses::bus] == 1;
        success *= result.bus_fault[0].disambiguation_string == "1";
        success *= result.bus_fault[0].monitored_variables.empty();

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
