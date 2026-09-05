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
  std::filesystem::remove_all(directory);

  Testing::TestingResults results;
  results += success.report("EMT study mu and constant signals");
  return results.summary();
}
