#include <iostream>

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
                "comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
                "freq_base": 60,
                "va_base": 100e6
            },
            "buses": [
                { "number": 1, "class": "bus", "name": "Bus 1", "init": {"Vr":0.994988, "Vi":0.099997}, "v_base": 115e3, "mon": ["Vr", "Vi"] },
                { "number": 2, "class": "infinite_bus", "name": "Bus 2", "init": {"Vr":1.0, "Vi":0}, "v_base": 115e3 }
            ],
            "devices": [
                { "class": "branch", "ports": {"bus1":1, "bus2":2}, "id": "1", "params": {"R":0, "X":0.1, "G":0, "B":0} },
                { "class": "GENROU", "ports": {"bus":1}, "id": "1", "params": {"p0":1, "q0":0.05013, "H":3, "D":0, "Ra":0, "Tdop":7, "Tdopp":0.04, "Tqopp":0.05,
                       "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqpp":0.18, "Xl":0.15, "S10":0, "S12":0}, "mon": ["delta", "omega"] }
                { "class": "bus_fault", "ports": {"bus":1}, "id": "1", "params": {"state0":0, "R":0, "X":1e-3} }
            ]
        })";

        TestStatus                success = true;
        GridDynamicsFormatHandler handler;
        Reader                    reader;
        StringStream              ss(json);

        reader.Parse(ss, handler);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
