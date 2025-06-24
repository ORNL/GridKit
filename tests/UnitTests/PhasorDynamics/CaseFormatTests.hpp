#include <iostream>

#include <rapidjson/error/en.h>

#include <Model/PhasorDynamics/GridDynamicsParser.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename RealT, typename IdxT>
    class CaseFormatTests
    {
      using FormatHandler = GridDynamicsFormatHandler<RealT, IdxT>;

    public:
      CaseFormatTests()  = default;
      ~CaseFormatTests() = default;

      TestOutcome simpleParse()
      {
        const char json[] =
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

        TestStatus    success = true;
        FormatHandler handler;
        Reader        reader;
        StringStream  ss(json);

        auto result  = reader.Parse(ss, handler);
        success     *= result;
        if (!result)
        {
          ParseErrorCode e = reader.GetParseErrorCode();
          size_t         o = reader.GetErrorOffset();
          std::cout << "parse error: " << GetParseError_En(e) << std::endl;
          std::cout << " at offset " << o << " near '" << std::string(json).substr(o, 20) << "...' with state " << handler.state << std::endl;
        }

        success *= handler.system_model.format_version == 0;
        success *= handler.system_model.format_revision == 1;
        success *= handler.system_model.case_name == "Two-bus test case 1";
        success *= handler.system_model.case_date_time == std::nullopt;
        success *= handler.system_model.case_description == "A two-bus test case for demonstrating the dynamics format";
        success *= handler.system_model.case_comments == "This case is set up to monitor the voltage at both buses and the machine angle and speed";
        success *= handler.system_model.freq_base == 60.0;
        success *= handler.system_model.va_base == 100.0e6;

        success *= handler.system_model.bus.size() == 2;
        success *= handler.system_model.branch.size() == 1;
        success *= handler.system_model.bus_fault.size() == 1;
        success *= handler.system_model.genrou.size() == 1;
        success *= handler.system_model.load.size() == 0;

        success *= handler.system_model.bus[0].bus_id == 1;
        success *= handler.system_model.bus[0].bus_type == FormatHandler::BusDataT::BusType::Default;
        success *= handler.system_model.bus[0].name == "Bus 1";
        success *= handler.system_model.bus[0].Vr0 == 0.994988;
        success *= handler.system_model.bus[0].Vi0 == 0.099997;
        success *= handler.system_model.bus[0].v_base == 115e3;
        success *= handler.system_model.bus[0].monitored_variables[static_cast<size_t>(FormatHandler::BusDataT::MonitorableVariables::Vr)];
        success *= handler.system_model.bus[0].monitored_variables[static_cast<size_t>(FormatHandler::BusDataT::MonitorableVariables::Vi)];
        success *= handler.system_model.bus[1].bus_id == 2;
        success *= handler.system_model.bus[1].bus_type == FormatHandler::BusDataT::BusType::Slack;
        success *= handler.system_model.bus[1].name == "Bus 2";
        success *= handler.system_model.bus[1].Vr0 == 1.0;
        success *= handler.system_model.bus[1].Vi0 == 0.0;
        success *= handler.system_model.bus[1].v_base == 115e3;
        success *= handler.system_model.bus[1].monitored_variables.none();

        success *= handler.system_model.branch[0].R == 0.0;
        success *= handler.system_model.branch[0].X == 0.1;
        success *= handler.system_model.branch[0].G == 0.0;
        success *= handler.system_model.branch[0].B == 0.0;
        success *= handler.system_model.branch[0].bus1_id == 1;
        success *= handler.system_model.branch[0].bus2_id == 2;
        success *= handler.system_model.branch[0].disambiguation_string == "1";
        success *= handler.system_model.branch[0].monitored_variables.none();

        success *= handler.system_model.genrou[0].unit_id == 1;
        success *= handler.system_model.genrou[0].p0 == 1.0;
        success *= handler.system_model.genrou[0].q0 == 0.05013;
        success *= handler.system_model.genrou[0].H == 3.0;
        success *= handler.system_model.genrou[0].D == 0.0;
        success *= handler.system_model.genrou[0].Ra == 0.0;
        success *= handler.system_model.genrou[0].Tdop == 7.0;
        success *= handler.system_model.genrou[0].Tdopp == 0.04;
        success *= handler.system_model.genrou[0].Tqop == 0.75;
        success *= handler.system_model.genrou[0].Tqopp == 0.05;
        success *= handler.system_model.genrou[0].Xd == 2.1;
        success *= handler.system_model.genrou[0].Xdp == 0.2;
        success *= handler.system_model.genrou[0].Xdpp == 0.18;
        success *= handler.system_model.genrou[0].Xq == 0.5;
        success *= handler.system_model.genrou[0].Xqp == 0.0;
        success *= handler.system_model.genrou[0].Xqpp == 0.18;
        success *= handler.system_model.genrou[0].Xl == 0.15;
        success *= handler.system_model.genrou[0].S10 == 0.0;
        success *= handler.system_model.genrou[0].S12 == 0.0;
        success *= handler.system_model.genrou[0].bus_id == 1;
        success *= handler.system_model.genrou[0].disambiguation_string == "1";
        success *= handler.system_model.genrou[0].monitored_variables[static_cast<size_t>(FormatHandler::GenrouDataT::MonitorableVariables::Delta)];
        success *= handler.system_model.genrou[0].monitored_variables[static_cast<size_t>(FormatHandler::GenrouDataT::MonitorableVariables::Omega)];

        success *= handler.system_model.bus_fault[0].R == 0.0;
        success *= handler.system_model.bus_fault[0].X == 1e-3;
        success *= handler.system_model.bus_fault[0].status == false;
        success *= handler.system_model.bus_fault[0].bus_id == 1;
        success *= handler.system_model.bus_fault[0].disambiguation_string == "1";
        success *= handler.system_model.bus_fault[0].monitored_variables.none();

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
