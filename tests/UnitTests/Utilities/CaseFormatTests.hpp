#include <iostream>

#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>
#include <nlohmann/json.hpp>

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
        using BusData = BusData<RealT, IdxT>;
        using BusType = BusData::BusType;

        const char data[] =
            R"({
               "header": {
                   "format_version": 0,
                   "format_revision": 1,
                   "case_name": "Two-bus test case 1",
                   "case_description": "A two-bus test case for demonstrating the dynamics format",
                   "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
                   "freq_base": 60.0,
                   "va_base": 100e6
               },
               "buses": [
                   { "number": 1, "class": "bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "v_base": 115e3, "mon": ["Vr", "Vi"] },
                   { "number": 2, "class": "infinite_bus", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0.0}, "v_base": 115e3 }
               ],
               "devices": [
                   { "class": "Branch", "ports": {"bus1":1, "bus2":2}, "id": "1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0} },
                   { "class": "Genrou", "ports": {"bus":1}, "id": "1", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05,
                          "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                   { "class": "bus_fault", "ports": {"bus":1}, "id": "1", "params": {"state0": false, "R":0.0, "X":1e-3} }
               ]
            })";

        TestStatus                   success = true;
        SystemModelData<RealT, IdxT> result  = json::parse(data);

        success *= result.format_version == 0;
        success *= result.format_revision == 1;
        success *= result.case_name == "Two-bus test case 1";
        success *= result.case_date_time == std::nullopt;
        success *= result.case_description == "A two-bus test case for demonstrating the dynamics format";
        success *= result.case_comments == "This case is set up to monitor the voltage at both buses and the machine angle and speed";
        success *= result.freq_base == 60.0;
        success *= result.va_base == 100.0e6;

        success *= result.bus.size() == 2;
        success *= result.branch.size() == 1;
        success *= result.bus_fault.size() == 1;
        success *= result.genrou.size() == 1;
        success *= result.load.size() == 0;

        success *= result.bus[0].bus_id == 1;
        success *= result.bus[0].bus_type == BusType::DEFAULT;
        success *= result.bus[0].name == "Bus 1";
        success *= result.bus[0].Vr0 == 0.994988;
        success *= result.bus[0].Vi0 == 0.099997;
        success *= result.bus[0].v_base == 115e3;
        success *= result.bus[0].monitored_variables[static_cast<size_t>(BusData::MonitorableVariables::VR)];
        success *= result.bus[0].monitored_variables[static_cast<size_t>(BusData::MonitorableVariables::VI)];
        success *= result.bus[1].bus_id == 2;
        success *= result.bus[1].bus_type == BusType::SLACK;
        success *= result.bus[1].name == "Bus 2";
        success *= result.bus[1].Vr0 == 1.0;
        success *= result.bus[1].Vi0 == 0.0;
        success *= result.bus[1].v_base == 115e3;
        success *= result.bus[1].monitored_variables.none();

        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::R]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::X]) == 0.1;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::G]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::B]) == 0.0;
        success *= result.branch[0].ports[BranchPorts::bus1] == 1;
        success *= result.branch[0].ports[BranchPorts::bus2] == 2;
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
        success *= result.genrou[0].ports[GenrouPorts::bus] == 1;
        success *= result.genrou[0].disambiguation_string == "1";
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::delta);
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::omega);

        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::R]) == 0.0;
        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::X]) == 1e-3;
        success *= !std::get<bool>(result.bus_fault[0].parameters[BusFaultParameters::state0]);
        success *= result.bus_fault[0].ports[BusFaultPorts::bus] == 1;
        success *= result.bus_fault[0].disambiguation_string == "1";
        success *= result.bus_fault[0].monitored_variables.empty();

        return success.report(__func__);
      }


      TestOutcome signalParse()
      {
        using namespace GridKit::PhasorDynamics;
        using BusData = BusData<RealT, IdxT>;
        using BusType = BusData::BusType;

        const char data[] =
            R"({
               "header": {
                   "format_version": 0,
                   "format_revision": 1,
                   "case_name": "Two-bus test case 1",
                   "case_description": "A two-bus test case for demonstrating the dynamics format",
                   "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
                   "freq_base": 60.0,
                   "va_base": 100e6
               },
               "buses": [
                   { "number": 1, "class": "bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "v_base": 115e3, "mon": ["Vr", "Vi"] },
                   { "number": 2, "class": "infinite_bus", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0.0}, "v_base": 115e3 }
               ],
               "signals": [
                   { "signal_id": 1, "name": "Machine Speed Deviation"},
                   { "signal_id": 2, "name": "Mechanical Power"},
                   { "signal_id": 3, "name": "Excitation Field"}
               ],
               "devices": [
                   { "class": "Branch", "ports": {"bus1":1, "bus2":2}, "id": "BR1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0} },
                   { "class": "Genrou", "ports": {"bus":1, "speed": 1, "pmech":2, "efd":3}, "id": "DV1", "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05, "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                   { "class": "Tgov1", "ports": {"bus":1, "speed": 1, "pmech":2}, "id": "DV2", "params": {"R":0.05, "T1":0.5,"T2":2.5, "T3":7.5, "Pvmax":0, "Pvmin":1, "Dt":0}},
                   { "class": "bus_fault", "ports": {"bus":1}, "id": "1", "params": {"state0": false, "R":0.0, "X":1e-3} }
               ]
            })";

        TestStatus                   success = true;
        SystemModelData<RealT, IdxT> result  = json::parse(data);

        success *= result.format_version == 0;
        success *= result.format_revision == 1;
        success *= result.case_name == "Two-bus test case 1";
        success *= result.case_date_time == std::nullopt;
        success *= result.case_description == "A two-bus test case for demonstrating the dynamics format";
        success *= result.case_comments == "This case is set up to monitor the voltage at both buses and the machine angle and speed";
        success *= result.freq_base == 60.0;
        success *= result.va_base == 100.0e6;

        success *= result.bus.size() == 2;
        success *= result.branch.size() == 1;
        success *= result.bus_fault.size() == 1;
        success *= result.genrou.size() == 1;
        success *= result.load.size() == 0;

        success *= result.bus[0].bus_id == 1;
        success *= result.bus[0].bus_type == BusType::DEFAULT;
        success *= result.bus[0].name == "Bus 1";
        success *= result.bus[0].Vr0 == 0.994988;
        success *= result.bus[0].Vi0 == 0.099997;
        success *= result.bus[0].v_base == 115e3;
        success *= result.bus[0].monitored_variables[static_cast<size_t>(BusData::MonitorableVariables::VR)];
        success *= result.bus[0].monitored_variables[static_cast<size_t>(BusData::MonitorableVariables::VI)];
        success *= result.bus[1].bus_id == 2;
        success *= result.bus[1].bus_type == BusType::SLACK;
        success *= result.bus[1].name == "Bus 2";
        success *= result.bus[1].Vr0 == 1.0;
        success *= result.bus[1].Vi0 == 0.0;
        success *= result.bus[1].v_base == 115e3;
        success *= result.bus[1].monitored_variables.none();

        success *= result.signal[0].signal_id == 1;
        success *= result.signal[0].name == "Machine Speed Deviation";
        success *= result.signal[1].signal_id == 2;
        success *= result.signal[1].name == "Mechanical Power";
        success *= result.signal[2].signal_id == 3;
        success *= result.signal[2].name == "Excitation Field";

        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::R]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::X]) == 0.1;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::G]) == 0.0;
        success *= std::get<RealT>(result.branch[0].parameters[BranchParameters::B]) == 0.0;
        success *= result.branch[0].ports[BranchPorts::bus1] == 1;
        success *= result.branch[0].ports[BranchPorts::bus2] == 2;
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
        success *= result.genrou[0].ports[GenrouPorts::bus] == 1;
        success *= result.genrou[0].disambiguation_string == "1";
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::delta);
        success *= result.genrou[0].monitored_variables.contains(GenrouMonitorableVariables::omega);

        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::R]) == 0.0;
        success *= std::get<RealT>(result.bus_fault[0].parameters[BusFaultParameters::X]) == 1e-3;
        success *= !std::get<bool>(result.bus_fault[0].parameters[BusFaultParameters::state0]);
        success *= result.bus_fault[0].ports[BusFaultPorts::bus] == 1;
        success *= result.bus_fault[0].disambiguation_string == "1";
        success *= result.bus_fault[0].monitored_variables.empty();

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
