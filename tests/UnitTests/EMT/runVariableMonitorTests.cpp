#include <cmath>
#include <limits>
#include <sstream>

#include <nlohmann/json.hpp>

#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  template <typename, typename>
  struct MonitorData
  {
    enum class MonitorableVariables
    {
      value,
      count,
      text
    };
  };

  template <typename, typename>
  struct Monitored
  {
  };

  using Monitor    = GridKit::Model::VariableMonitor<Monitored<double, size_t>, MonitorData>;
  using Controller = GridKit::Model::VariableMonitorController<double>;
  using Variable   = MonitorData<double, size_t>::MonitorableVariables;
  using Format     = GridKit::Model::VariableMonitorFormat;
  using GridKit::Testing::TestStatus;

  bool formatting()
  {
    bool success = true;
    for (double value : {-0.0, 0.0, 1.25, -1.25, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      std::string actual = "prefix:";
      GridKit::Model::VariableMonitorDetail::appendReal(actual, value);
      std::ostringstream expected;
      expected.precision(std::numeric_limits<double>::digits10 + 1);
      expected << "prefix:" << std::scientific << value;
      success &= actual == expected.str();
    }
    double  time = 0.0;
    Monitor monitor("sample", {Variable::count, Variable::text});
    monitor.set(Variable::count, []
                { return 17; });
    monitor.set(Variable::text, []
                { return "ready"; });
    Controller controller(time);
    controller.addMonitor(&monitor);
    std::ostringstream output;
    controller.printFull(output, Controller::Csv{});
    success &= output.str() == "0.0000000000000000e+00,17,ready\n";
    return success;
  }

  bool sinks()
  {
    double  time = 0.0, value = 1.25;
    int     reads = 0;
    Monitor monitor("sample", {Variable::value, Variable::count});
    monitor.set(Variable::value, [&]
                { ++reads; return value; });
    monitor.set(Variable::count, []
                { return 17; });
    Controller controller(time);
    controller.addMonitor(&monitor);
    std::ostringstream csv1, csv2, semicolon, json1, json2, yaml1, yaml2;
    controller.addSink({Format::CSV}, csv1);
    controller.addSink({Format::JSON}, json1);
    controller.addSink({Format::CSV, {}, ";"}, semicolon);
    controller.addSink({Format::CSV}, csv2);
    controller.addSink({Format::JSON}, json2);
    controller.addSink({Format::YAML}, yaml1);
    controller.addSink({Format::YAML}, yaml2);
    controller.start();
    controller.print();
    // Same timestamp, new state: every format must sample again.
    value = -2.5;
    controller.print();
    controller.stop();
    bool success  = reads == 8;
    success      &= csv1.str() == csv2.str();
    success      &= json1.str() == json2.str();
    success      &= yaml1.str() == yaml2.str();
    auto comma    = csv1.str();
    std::replace(comma.begin(), comma.end(), ',', ';');
    success         &= comma == semicolon.str();
    const auto json  = nlohmann::json::parse(json1.str());
    success         &= json.size() == 2;
    success         &= json[0]["sample"]["value"] == 1.25;
    success         &= json[1]["sample"]["value"] == -2.5;
    success         &= json[0]["sample"]["count"] == 17;
    return success;
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  TestStatus                       format  = formatting();
  results                                 += format.report("monitor scalar formatting");
  TestStatus output                        = sinks();
  results                                 += output.report("monitor sink reuse and row lifecycle");
  return results.summary();
}
