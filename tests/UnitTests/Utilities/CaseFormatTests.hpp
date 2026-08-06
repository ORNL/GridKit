#pragma once

#include <iostream>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
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
        using ReecbData = Controller::ReecbData<RealT, IdxT>;

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
                   { "class": "Regca", "ports": {"bus":1}, "id": "CV1", "params": {"p0":0.0, "q0":0.0, "mva":100, "Tg":0.02, "TM":0.02, "Rqmax":999.0, "Rqmin":-999.0, "Rpmax":999.0, "sL":true, "IL1":1.1, "VL0":0.4, "VL1":0.9, "VA0":0.4, "VA1":0.9, "Vhvmax":1.2}, "mon": ["ir", "ii", "p", "q"] },
                   { "class": "Reecb", "ports": {"bus":1}, "id": "REE1", "params": {"mva":50.0, "Pqflag":true}, "mon": ["iqcmd", "ipcmd", "vmeas", "pmeas"] },
                   { "class": "BusFault", "ports": {"bus":1}, "id": "1", "params": {"state0": false, "R":0.0, "X":1e-3} }
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
        success *= result.reecb.size() == 1;
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
        success *= std::get<RealT>(result.reecb[0].parameters[ReecbData::Parameters::mva]) == 50.0;
        success *= std::get<bool>(result.reecb[0].parameters[ReecbData::Parameters::Pqflag]);
        success *= result.reecb[0].buses[ReecbData::Buses::bus] == 1;
        success *= result.reecb[0].disambiguation_string == "REE1";
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::iqcmd);
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::ipcmd);
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::vmeas);
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::pmeas);

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
        using BusData     = BusData<RealT, IdxT>;
        using BusType     = typename BusData::BusType;
        using Esdc1aData  = Exciter::Esdc1aData<RealT, IdxT>;
        using GastPtiData = Governor::GastPtiData<RealT, IdxT>;
        using HygovData   = Governor::HygovData<RealT, IdxT>;
        using RepcaData   = Controller::RepcaData<RealT, IdxT>;

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
                   { "signal_id": 6, "name": "Under-excitation Limiter"},
                   { "signal_id": 7, "name": "Hydro Mechanical Power"},
                   { "signal_id": 8, "name": "Governor Load Reference"},
                   { "signal_id": 9, "name": "Governor Auxiliary Power"},
                   { "signal_id": 10, "name": "Governor Reference"},
                   { "signal_id": 11, "name": "Branch Current Real"},
                   { "signal_id": 12, "name": "Branch Current Imaginary"},
                   { "signal_id": 13, "name": "Branch Active Power"},
                   { "signal_id": 14, "name": "Branch Reactive Power"},
                   { "signal_id": 15, "name": "Frequency"},
                   { "signal_id": 16, "name": "Plant Voltage Reference"},
                   { "signal_id": 17, "name": "Plant Active Power Reference"},
                   { "signal_id": 18, "name": "Reactive Power Reference"},
                   { "signal_id": 19, "name": "Frequency Reference"},
                   { "signal_id": 20, "name": "Reactive Power Command"},
                   { "signal_id": 21, "name": "Active Power Command"},
                   { "signal_id": 22, "name": "Electrical Power"},
                   { "signal_id": 23, "name": "Reactive Power"},
                   { "signal_id": 24, "name": "Reactive Reference"},
                   { "signal_id": 25, "name": "Power Factor Reference"},
                   { "signal_id": 26, "name": "Active Power Reference"},
                   { "signal_id": 27, "name": "Reactive Current Command"},
                   { "signal_id": 28, "name": "Active Current Command"}
               ],
               "devices": [
                   { "class": "Branch", "ports": {"bus1":1, "bus2":2}, "id": "BR1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0, "tap":1.05, "phase":0.1} },
                   { "class": "Genrou", "ports": {"bus":1, "speed": 1, "pmech":2, "efd":3}, "id": "DV1", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05, "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                   { "class": "Tgov1", "ports": {"speed": 1, "pref":10, "pmech":2}, "id": "DV2", "params": {"R":0.05, "T1":0.5,"T2":2.5, "T3":7.5, "Pvmax":0.0, "Pvmin":1.0, "Dt":0.0}},
                   { "class": "Esdc1a", "ports": {"bus":1, "speed":1, "vref":4, "vs":5, "vuel":6, "efd":3}, "id": "DV5", "params": {"Tr":0.0, "Ka":40.0, "Ta":0.1, "Tb":0.0, "Tc":0.0, "Vrmax":1.0, "Vrmin":-1.0, "Ke":0.1, "Te":0.5, "Kf":0.05, "Tf1":0.7, "Spdmlt":false, "E1":2.8, "Se1":0.08, "E2":3.7, "Se2":0.33, "UEL":0, "exclim":true}, "mon": ["efd", "vc", "vr", "vf", "se", "vfe"] },
                   { "class": "Hygov", "ports": {"speed": 1, "pmech": 7, "pref": 8, "paux": 9}, "id": "DV6", "params": {"Trate": 80.0, "Rperm": 0.05, "Rtemp": 0.35, "Tr": 5.0, "Tf": 0.05, "Tg": 0.5,
                          "Velm": 0.2, "Gmax": 0.98, "Gmin": 0.02, "Tw": 1.2, "At": 1.1, "Dturb": 0.4, "Qnl": 0.08, "Tn": 0.7, "Tnp": 1.4, "db1": 0.01, "db2": 0.02, "Hdam": 1.05,
                          "Gv0": 0.0, "Gv1": 0.2, "Gv2": 0.4, "Gv3": 0.6, "Gv4": 0.8, "Gv5": 1.0,
                          "Pgv0": 0.0, "Pgv1": 0.15, "Pgv2": 0.42, "Pgv3": 0.66, "Pgv4": 0.85, "Pgv5": 1.0}, "mon": ["pmech", "filter", "desiredgate", "gate", "flow", "head"]},
                   { "class": "Repca", "ports": {"bus":1, "ir":11, "ii":12, "p":13, "q":14, "freq":15, "vref":16, "pref":17, "qref":18, "freqref":19, "qext":20, "pext":21}, "id": "PC1", "params": {"mva":50, "VcompFlag":false, "RefFlag":true, "Freqflag":true, "Tfltr":0.2, "Vfrz":0.65, "Rc":0.02, "Xc":0.03, "Kc":0.4, "dbdlow":-0.02, "dbdupper":0.03, "emax":0.8, "emin":-0.7, "Kp":2.0, "Ki":3.0, "Qmax":0.9, "Qmin":-0.8, "Tft":0.2, "Tfv":1.5, "Tp":0.4, "fdbd1":-0.01, "fdbd2":0.015, "Ddn":2.0, "Dup":1.0, "femax":0.6, "femin":-0.5, "Kpg":1.7, "Kig":1.8, "Pmax":1.2, "Pmin":0.1, "Tlag":0.5}, "mon": ["qext", "pext", "vmeas", "qmeas", "pmeas"] },
                   { "class": "Ieeet1", "ports": {"bus":1, "speed": 1, "efd":3}, "id": "DV3", "params": {"Tr":0.0, "Ka":50.0, "Ta":0.04, "Ke":-0.06, "Te":0.6, "Kf":0.09, "Tf":1.46, "Vrmin":-1.0, "Vrmax":1.0, "E1":2.8, "E2":3.373, "Se1":0.04, "Se2":0.33, "Ispdlim":0.0}},
                   { "class": "SexsPti", "ports": {"bus":1, "efd":3}, "id": "DV4", "params": {"Ta":0.1, "Tb":0.5, "Te":0.8, "K":10.0, "Efdmax":5.0, "Efdmin":-5.0}},
                   { "class": "GastPti", "ports": {"speed":1, "pmech":2, "pref":10}, "id": "DV7", "params": {"R":0.045, "T1":0.42, "T2":0.12, "T3":3.2, "At":0.95, "Kt":2.2, "Vmax":1.05, "Vmin":0.15, "Dturb":0.02, "Trate":120.0}, "mon": ["pmech", "xvalve", "xflow", "xtemp", "vload", "vtemp"] },
                   { "class": "Reecb", "ports": {"bus":1, "pe":22, "qgen":23, "qext":24, "pfaref":25, "pref":26, "iqcmd":27, "ipcmd":28}, "id": "EC1", "params": {"mva":100.0, "PfFlag":false, "VFlag":true, "QFlag":true, "Pqflag":true, "Trv":0.02, "Tp":0.05, "Vref0":1.0, "Vdip":0.85, "Vup":1.15, "dbd1":-0.01, "dbd2":0.01, "kqv":5.0, "Iql1":-1.1, "Iqh1":1.1, "Qmax":0.436, "Qmin":-0.436, "Kqp":0.1, "Kqi":0.2, "Vmax":1.1, "Vmin":0.9, "Kvp":18.0, "Kvi":5.0, "Tiq":0.02, "Tpord":0.02, "dPmax":99.0, "dPmin":-99.0, "Pmax":1.0, "Pmin":0.0, "Imax":1.3}, "mon": ["iqcmd", "ipcmd", "vmeas", "pmeas"]},
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
        success *= result.gastpti.size() == 1;
        success *= result.hygov.size() == 1;
        success *= result.repca.size() == 1;
        success *= result.loadz.size() == 0;
        success *= result.exciter.size() == 1;
        success *= result.sexspti.size() == 1;
        success *= result.reecb.size() == 1;
        success *= result.signal.size() == 28;

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
        success *= result.signal[6].signal_id == 7;
        success *= result.signal[6].name == "Hydro Mechanical Power";
        success *= result.signal[7].signal_id == 8;
        success *= result.signal[7].name == "Governor Load Reference";
        success *= result.signal[8].signal_id == 9;
        success *= result.signal[8].name == "Governor Auxiliary Power";
        success *= result.signal[9].signal_id == 10;
        success *= result.signal[9].name == "Governor Reference";
        success *= result.signal[27].signal_id == 28;
        success *= result.signal[27].name == "Active Current Command";

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
        success *= result.gov[0].signal_inputs[Governor::Tgov1SignalInputs::pref] == 10;
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

        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::R]) == 0.045;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::T1]) == 0.42;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::T2]) == 0.12;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::T3]) == 3.2;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::At]) == 0.95;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::Kt]) == 2.2;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::Vmax]) == 1.05;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::Vmin]) == 0.15;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::Dturb]) == 0.02;
        success *= std::get<RealT>(result.gastpti[0].parameters[GastPtiData::Parameters::Trate]) == 120.0;
        success *= result.gastpti[0].signal_inputs[GastPtiData::SignalInputs::speed] == 1;
        success *= result.gastpti[0].signal_inputs[GastPtiData::SignalInputs::pref] == 10;
        success *= result.gastpti[0].signal_outputs[GastPtiData::SignalOutputs::pmech] == 2;
        success *= result.gastpti[0].disambiguation_string == "DV7";
        success *= result.gastpti[0].monitored_variables.contains(GastPtiData::MonitorableVariables::pmech);
        success *= result.gastpti[0].monitored_variables.contains(GastPtiData::MonitorableVariables::xvalve);
        success *= result.gastpti[0].monitored_variables.contains(GastPtiData::MonitorableVariables::xflow);
        success *= result.gastpti[0].monitored_variables.contains(GastPtiData::MonitorableVariables::xtemp);
        success *= result.gastpti[0].monitored_variables.contains(GastPtiData::MonitorableVariables::vload);
        success *= result.gastpti[0].monitored_variables.contains(GastPtiData::MonitorableVariables::vtemp);

        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Trate]) == 80.0;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Rperm]) == 0.05;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Rtemp]) == 0.35;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Tr]) == 5.0;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Tf]) == 0.05;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Tg]) == 0.5;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Velm]) == 0.2;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gmax]) == 0.98;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gmin]) == 0.02;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Tw]) == 1.2;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::At]) == 1.1;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Dturb]) == 0.4;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Qnl]) == 0.08;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Tn]) == 0.7;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Tnp]) == 1.4;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::db1]) == 0.01;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::db2]) == 0.02;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Hdam]) == 1.05;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gv0]) == 0.0;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gv1]) == 0.2;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gv2]) == 0.4;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gv3]) == 0.6;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gv4]) == 0.8;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Gv5]) == 1.0;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Pgv0]) == 0.0;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Pgv1]) == 0.15;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Pgv2]) == 0.42;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Pgv3]) == 0.66;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Pgv4]) == 0.85;
        success *= std::get<RealT>(result.hygov[0].parameters[HygovData::Parameters::Pgv5]) == 1.0;
        success *= result.hygov[0].signal_inputs[HygovData::SignalInputs::speed] == 1;
        success *= result.hygov[0].signal_inputs[HygovData::SignalInputs::pref] == 8;
        success *= result.hygov[0].signal_inputs[HygovData::SignalInputs::paux] == 9;
        success *= result.hygov[0].signal_outputs[HygovData::SignalOutputs::pmech] == 7;
        success *= result.hygov[0].disambiguation_string == "DV6";
        success *= result.hygov[0].monitored_variables.contains(HygovData::MonitorableVariables::pmech);
        success *= result.hygov[0].monitored_variables.contains(HygovData::MonitorableVariables::filter);
        success *= result.hygov[0].monitored_variables.contains(
            HygovData::MonitorableVariables::desiredgate);
        success *= result.hygov[0].monitored_variables.contains(HygovData::MonitorableVariables::gate);
        success *= result.hygov[0].monitored_variables.contains(HygovData::MonitorableVariables::flow);
        success *= result.hygov[0].monitored_variables.contains(HygovData::MonitorableVariables::head);

        success *= std::get<IdxT>(result.repca[0].parameters[RepcaData::Parameters::mva]) == 50;
        success *= !std::get<bool>(result.repca[0].parameters[RepcaData::Parameters::VcompFlag]);
        success *= std::get<bool>(result.repca[0].parameters[RepcaData::Parameters::RefFlag]);
        success *= std::get<bool>(result.repca[0].parameters[RepcaData::Parameters::Freqflag]);
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Tfltr]) == 0.2;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Vfrz]) == 0.65;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Rc]) == 0.02;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Xc]) == 0.03;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Kc]) == 0.4;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::dbdlow]) == -0.02;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::dbdupper]) == 0.03;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::emax]) == 0.8;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::emin]) == -0.7;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Kp]) == 2.0;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Ki]) == 3.0;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Qmax]) == 0.9;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Qmin]) == -0.8;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Tft]) == 0.2;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Tfv]) == 1.5;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Tp]) == 0.4;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::fdbd1]) == -0.01;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::fdbd2]) == 0.015;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Ddn]) == 2.0;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Dup]) == 1.0;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::femax]) == 0.6;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::femin]) == -0.5;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Kpg]) == 1.7;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Kig]) == 1.8;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Pmax]) == 1.2;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Pmin]) == 0.1;
        success *= std::get<RealT>(result.repca[0].parameters[RepcaData::Parameters::Tlag]) == 0.5;
        success *= result.repca[0].buses[RepcaData::Buses::bus] == 1;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::ir] == 11;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::ii] == 12;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::p] == 13;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::q] == 14;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::freq] == 15;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::vref] == 16;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::pref] == 17;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::qref] == 18;
        success *= result.repca[0].signal_inputs[RepcaData::SignalInputs::freqref] == 19;
        success *= result.repca[0].signal_outputs[RepcaData::SignalOutputs::qext] == 20;
        success *= result.repca[0].signal_outputs[RepcaData::SignalOutputs::pext] == 21;
        success *= result.repca[0].disambiguation_string == "PC1";
        success *= result.repca[0].monitored_variables.contains(RepcaData::MonitorableVariables::qext);
        success *= result.repca[0].monitored_variables.contains(RepcaData::MonitorableVariables::pext);
        success *= result.repca[0].monitored_variables.contains(RepcaData::MonitorableVariables::vmeas);
        success *= result.repca[0].monitored_variables.contains(RepcaData::MonitorableVariables::qmeas);
        success *= result.repca[0].monitored_variables.contains(RepcaData::MonitorableVariables::pmeas);

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

        using ReecbData    = Controller::ReecbData<RealT, IdxT>;
        using ReecbParams  = ReecbData::Parameters;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::mva]) == 100.0;
        success           *= !std::get<bool>(result.reecb[0].parameters[ReecbParams::PfFlag]);
        success           *= std::get<bool>(result.reecb[0].parameters[ReecbParams::VFlag]);
        success           *= std::get<bool>(result.reecb[0].parameters[ReecbParams::QFlag]);
        success           *= std::get<bool>(result.reecb[0].parameters[ReecbParams::Pqflag]);
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Trv]) == 0.02;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Tp]) == 0.05;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Vref0]) == 1.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Vdip]) == 0.85;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Vup]) == 1.15;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::dbd1]) == -0.01;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::dbd2]) == 0.01;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::kqv]) == 5.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Iql1]) == -1.1;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Iqh1]) == 1.1;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Qmax]) == 0.436;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Qmin]) == -0.436;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Kqp]) == 0.1;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Kqi]) == 0.2;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Vmax]) == 1.1;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Vmin]) == 0.9;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Kvp]) == 18.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Kvi]) == 5.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Tiq]) == 0.02;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Tpord]) == 0.02;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::dPmax]) == 99.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::dPmin]) == -99.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Pmax]) == 1.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Pmin]) == 0.0;
        success           *= std::get<RealT>(result.reecb[0].parameters[ReecbParams::Imax]) == 1.3;
        success           *= result.reecb[0].buses[ReecbData::Buses::bus] == 1;
        success           *= result.reecb[0].signal_inputs[ReecbData::SignalInputs::pe] == 22;
        success           *= result.reecb[0].signal_inputs[ReecbData::SignalInputs::qgen] == 23;
        success           *= result.reecb[0].signal_inputs[ReecbData::SignalInputs::qext] == 24;
        success           *= result.reecb[0].signal_inputs[ReecbData::SignalInputs::pfaref] == 25;
        success           *= result.reecb[0].signal_inputs[ReecbData::SignalInputs::pref] == 26;
        success           *= result.reecb[0].signal_outputs[ReecbData::SignalOutputs::iqcmd] == 27;
        success           *= result.reecb[0].signal_outputs[ReecbData::SignalOutputs::ipcmd] == 28;
        success           *= result.reecb[0].disambiguation_string == "EC1";
        success           *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::iqcmd);
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::ipcmd);
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::vmeas);
        success *= result.reecb[0].monitored_variables.contains(
            ReecbData::MonitorableVariables::pmeas);

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
