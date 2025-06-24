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
    public:
      CaseFormatTests()  = default;
      ~CaseFormatTests() = default;

      TestOutcome simpleParse()
      {
        const char json[] = R"(
        {
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

        TestStatus                success = true;
        GridDynamicsFormatHandler handler;
        Reader                    reader;
        StringStream              ss(json);

        auto result  = reader.Parse(ss, handler);
        success     *= result;
        if (!result)
        {
          ParseErrorCode e = reader.GetParseErrorCode();
          size_t         o = reader.GetErrorOffset();
          std::cout << "parse error: " << GetParseError_En(e) << std::endl;
          std::cout << " at offset " << o << " near '" << std::string(json).substr(o, 20) << "...' with state " << handler.state << std::endl;
        }

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
