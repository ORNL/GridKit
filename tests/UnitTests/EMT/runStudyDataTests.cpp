#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>

#include <GridKit/Testing/Testing.hpp>

#include <application/EMT/AnalysisUtilities.hpp>

int main()
{
  using namespace GridKit;
  using json                  = nlohmann::json;
  Testing::TestStatus success = true;
  auto                rejects = [](auto&& action)
  {
    try
    {
      action();
    }
    catch (const std::exception&)
    {
      return true;
    }
    return false;
  };
  const json base      = {{"system_model_file", "case.json"}, {"tmax", 1.0}};
  const auto defaults  = base.get<EMT::StudyData>();
  success             *= defaults.mu == Math::DEFAULT_MU<double>;
  auto input           = base;
  input["mu"]          = 50000.0;
  const auto sharp     = input.get<EMT::StudyData>();
  EMT::configureCommonMath<double>(sharp);
  success *= std::abs(Math::ramp(0.0) - std::log(2.0) / 50000.0) < 1e-18;
  EMT::configureCommonMath<double>(defaults);
  success *= Math::MU<double> == Math::DEFAULT_MU<double>;
  for (const auto& invalid : {json(0), json(-1), json(nullptr), json("50000"), json(true), json(std::numeric_limits<double>::infinity())})
  {
    input["mu"]  = invalid;
    success     *= rejects([&]
                       { input.get<EMT::StudyData>(); });
  }

  const json switch_event  = {{"time", 0.0}, {"type", "switch"}, {"element_id", "plant.breaker"}, {"open", false}};
  const json signal_event  = {{"time", 0.5}, {"type", "signal_step"}, {"signal_id", "plant.reference"}, {"value", 0.8}};
  input                    = base;
  input["events"]          = {switch_event, signal_event, signal_event};
  const auto events        = input.get<EMT::StudyData>().events;
  success                 *= events.size() == 3 && events[0].time == 0.0 && events[1].time == 0.5;
  success                 *= std::get<EMT::SwitchEvent>(events[0].action).element_id == "plant.breaker";
  success                 *= !std::get<EMT::SwitchEvent>(events[0].action).open;
  success                 *= std::get<EMT::SignalStep>(events[1].action).signal_id == "plant.reference";
  success                 *= std::get<EMT::SignalStep>(events[1].action).value == 0.8;
  auto rejects_event       = [&](const json& event)
  {
    input["events"] = {event};
    return rejects([&]
                   { input.get<EMT::StudyData>(); });
  };
  for (const auto& invalid : {json(true), json("0.5"), json(-0.1), json(1.1), json(std::numeric_limits<double>::infinity())})
  {
    auto event     = signal_event;
    event["time"]  = invalid;
    success       *= rejects_event(event);
  }
  for (const auto& invalid : {json(true), json("0.8"), json(nullptr), json(std::numeric_limits<double>::infinity())})
  {
    auto event      = signal_event;
    event["value"]  = invalid;
    success        *= rejects_event(event);
  }
  for (const auto& invalid : {json(0), json("false"), json(nullptr)})
  {
    auto event     = switch_event;
    event["open"]  = invalid;
    success       *= rejects_event(event);
  }
  for (const auto& type : {"switch_open", "switch_close", "unknown"})
  {
    auto event     = switch_event;
    event["type"]  = type;
    success       *= rejects_event(event);
  }
  for (const auto& key : {"time", "type", "signal_id", "value"})
  {
    auto event = signal_event;
    event.erase(key);
    success *= rejects_event(event);
  }
  auto extra           = signal_event;
  extra["element_id"]  = "plant.breaker";
  success             *= rejects_event(extra);
  input["events"]      = {signal_event, switch_event};
  success             *= rejects([&]
                     { input.get<EMT::StudyData>(); });
  input                = base;
  for (double invalid : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
  {
    input["tmax"]  = invalid;
    success       *= rejects([&]
                       { input.get<EMT::StudyData>(); });
  }

  const auto directory = std::filesystem::temp_directory_path()
                         / ("gridkit-emt-study-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(directory);
  const json model = {
      {"header", {{"case_name", "constants"}, {"case_description", "constant overrides"}, {"case_comments", ""}}},
      {"signals", {{{"id", "dc"}, {"value", 1000.0}}, {{"id", "output"}}}},
      {"devices", {{{"class", "Container"}, {"id", "child"}, {"signals", {{{"id", "dc"}, {"value", 2000.0}}}}, {"devices", json::array()}}}}};
  std::ofstream(directory / "case.json") << model.dump();
  input                  = base;
  input["signal_values"] = {{"dc", 3000.0}, {"child.dc", 4000.0}};
  auto parse             = [&]
  {
    std::ofstream(directory / "study.json") << input.dump();
    return EMT::parseStudyData(directory / "study.json");
  };
  const auto overridden  = parse();
  success               *= overridden.model_data.signal[0].value == 3000.0;
  success               *= overridden.model_data.container[0].signal[0].value == 4000.0;
  for (const auto& name : {"output", "missing", "missing.dc", "child.missing"})
  {
    input["signal_values"]  = {{name, 5.0}};
    success                *= rejects(parse);
  }
  input["signal_values"]  = json::array({3.0});
  success                *= rejects(parse);
  input["signal_values"]  = {{"dc", "3000"}};
  success                *= rejects(parse);

  input               = base;
  input["state_file"] = "state.json";
  auto parse_state    = [&](const json& state)
  {
    std::ofstream(directory / "state.json") << state.dump();
    return parse().state;
  };
  const auto state  = parse_state({{"buses", {{"child.bus", {{"va", 10}, {"vb", -5.0}, {"vc", nullptr}, {"injections", json::array()}}}}},
                                   {"devices", {{"child.load", {{"ia", -2.0}}}, {"switch", {{"open", true}}}, {"unused", nullptr}}}});
  success          *= state == std::map<std::string, std::map<std::string, double>>{{"child.bus", {{"va", 10.0}, {"vb", -5.0}}}, {"child.load", {{"ia", -2.0}}}, {"switch", {{"open", 1.0}}}};
  success          *= parse_state({{"devices", {{"load", {{"ia", nullptr}}}}}}).empty();
  for (const auto& invalid : {json(true), json("1"), json::array({1})})
  {
    success *= rejects([&]
                       { parse_state({{"devices", {{"load", {{"ia", invalid}}}}}}); });
  }
  for (const auto& invalid : {json::array(), json{{"buses", 1}}, json{{"devices", {{"load", 1}}}}})
  {
    success *= rejects([&]
                       { parse_state(invalid); });
  }
  std::filesystem::remove_all(directory);

  Testing::TestingResults results;
  results += success.report("EMT study mu, constants, typed events, and state");
  return results.summary();
}
