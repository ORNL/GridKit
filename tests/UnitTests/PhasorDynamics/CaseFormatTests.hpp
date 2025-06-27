#include <iostream>

#include <Model/PhasorDynamics/SystemModelData.hpp>
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
      using SystemModelDataT = PhasorDynamics::SystemModelData<RealT, IdxT>;

    public:
      CaseFormatTests()  = default;
      ~CaseFormatTests() = default;

      TestOutcome simpleParse()
      {
        const char data[] =
            R"({
                   "header": {
                       "format_version": 0,
                       "format_revision": 1,
                       "case_name": "Two-bus test case 1",
                       "case_description": "A two-bus test case for demonstrating the dynamics format",
                       "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
                       "freq_base": 60.0,
                       "va_base": 100.0e6
                   },
                   "buses": [
                       { "number": 1, "class": "bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "v_base": 115e3, "mon": ["Vr", "Vi"] },
                       { "number": 2, "class": "infinite_bus", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0.0}, "v_base": 115e3 }
                   ],
                   "devices": [
                       { "class": "branch", "ports": {"bus1":1, "bus2":2}, "id": "1", "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0} },
                       { "class": "GENROU", "ports": {"bus":1}, "id": "1", "params": {"unit_id": 1, "p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05,
                              "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.0, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0}, "mon": ["delta", "omega"] },
                       { "class": "bus_fault", "ports": {"bus":1}, "id": "1", "params": {"state0": false, "R":0.0, "X":1e-3} }
                   ]
               })";

        TestStatus       success = true;
        SystemModelDataT result  = json::parse(data);

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
        success *= result.bus[0].bus_type == SystemModelDataT::BusDataT::BusType::DEFAULT;
        success *= result.bus[0].name == "Bus 1";
        success *= result.bus[0].Vr0 == 0.994988;
        success *= result.bus[0].Vi0 == 0.099997;
        success *= result.bus[0].v_base == 115e3;
        success *= result.bus[0].monitored_variables[static_cast<size_t>(SystemModelDataT::BusDataT::MonitorableVariables::VR)];
        success *= result.bus[0].monitored_variables[static_cast<size_t>(SystemModelDataT::BusDataT::MonitorableVariables::VI)];
        success *= result.bus[1].bus_id == 2;
        success *= result.bus[1].bus_type == SystemModelDataT::BusDataT::BusType::SLACK;
        success *= result.bus[1].name == "Bus 2";
        success *= result.bus[1].Vr0 == 1.0;
        success *= result.bus[1].Vi0 == 0.0;
        success *= result.bus[1].v_base == 115e3;
        success *= result.bus[1].monitored_variables.none();

        success *= result.branch[0].R == 0.0;
        success *= result.branch[0].X == 0.1;
        success *= result.branch[0].G == 0.0;
        success *= result.branch[0].B == 0.0;
        success *= result.branch[0].bus1_id == 1;
        success *= result.branch[0].bus2_id == 2;
        success *= result.branch[0].disambiguation_string == "1";
        success *= result.branch[0].monitored_variables.none();

        success *= result.genrou[0].unit_id == 1;
        success *= result.genrou[0].p0 == 1.0;
        success *= result.genrou[0].q0 == 0.05013;
        success *= result.genrou[0].H == 3.0;
        success *= result.genrou[0].D == 0.0;
        success *= result.genrou[0].Ra == 0.0;
        success *= result.genrou[0].Tdop == 7.0;
        success *= result.genrou[0].Tdopp == 0.04;
        success *= result.genrou[0].Tqop == 0.75;
        success *= result.genrou[0].Tqopp == 0.05;
        success *= result.genrou[0].Xd == 2.1;
        success *= result.genrou[0].Xdp == 0.2;
        success *= result.genrou[0].Xdpp == 0.18;
        success *= result.genrou[0].Xq == 0.5;
        success *= result.genrou[0].Xqp == 0.0;
        success *= result.genrou[0].Xqpp == 0.18;
        success *= result.genrou[0].Xl == 0.15;
        success *= result.genrou[0].S10 == 0.0;
        success *= result.genrou[0].S12 == 0.0;
        success *= result.genrou[0].bus_id == 1;
        success *= result.genrou[0].disambiguation_string == "1";
        success *= result.genrou[0].monitored_variables[static_cast<size_t>(SystemModelDataT::GenrouDataT::MonitorableVariables::DELTA)];
        success *= result.genrou[0].monitored_variables[static_cast<size_t>(SystemModelDataT::GenrouDataT::MonitorableVariables::OMEGA)];

        success *= result.bus_fault[0].R == 0.0;
        success *= result.bus_fault[0].X == 1e-3;
        success *= result.bus_fault[0].status == false;
        success *= result.bus_fault[0].bus_id == 1;
        success *= result.bus_fault[0].disambiguation_string == "1";
        success *= result.bus_fault[0].monitored_variables.none();

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
