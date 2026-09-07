#include <cstdint>
#include <limits>
#include <string>

#include <GridKit/Model/EMT/Operators/Rational/StateSpace/StateSpaceDataJSONParser.hpp>
#include <GridKit/Model/EMT/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using json       = nlohmann::json;
  using SourceData = EMT::VoltageSourceData<double, size_t>;
  using Parameter  = SourceData::Parameters;

  template <typename Action>
  bool rejects(Action action, const std::string& context = {})
  {
    try
    {
      action();
    }
    catch (const std::exception& error)
    {
      return std::string(error.what()).find(context) != std::string::npos;
    }
    return false;
  }

  Testing::TestOutcome numericParameters()
  {
    Testing::TestStatus success = true;
    for (const auto& scalar : {json(3), json(3.0)})
    {
      const auto data  = json{{"class", "VoltageSource"}, {"id", "source"}, {"params", {{"omega", scalar}, {"N", scalar}, {"E", {1, 2, 3}}, {"phi", {-1, 0, 1}}}}}.get<SourceData>();
      success         *= EMT::parameter<double>(data, Parameter::omega) == 3.0;
      success         *= EMT::parameter<size_t>(data, Parameter::N) == 3;
      success         *= EMT::parameter<EMT::ABCVector<double>>(data, Parameter::E) == EMT::ABCVector<double>{1, 2, 3};
      success         *= EMT::parameter<EMT::ABCVector<double>>(data, Parameter::phi) == EMT::ABCVector<double>{-1, 0, 1};
      success         *= EMT::parameter<double>(data, Parameter::Rs, 0.25) == 0.25;
      success         *= rejects([&]
                         { EMT::parameter<double>(data, Parameter::Ls); },
                         "Ls");
    }
    const auto negative  = json{{"class", "VoltageSource"}, {"id", "negative"}, {"params", {{"omega", -2}}}}.get<SourceData>();
    success             *= EMT::parameter<double>(negative, Parameter::omega) == -2.0;
    success             *= rejects([&]
                       { EMT::parameter<size_t>(negative, Parameter::omega); },
                       "negative");

    using NarrowData = EMT::VoltageSourceData<double, uint8_t>;
    for (const auto& scalar : {json(0), json(0.0), json(255), json(255.0)})
    {
      const auto data  = json{{"class", "VoltageSource"}, {"id", "narrow"}, {"params", {{"N", scalar}}}}.get<NarrowData>();
      success         *= EMT::parameter<uint8_t>(data, Parameter::N) == scalar.get<uint8_t>();
    }
    for (const auto& scalar : {json(-1), json(-1.0), json(0.5), json(256), json(256.0), json(true)})
    {
      const auto data  = json{{"class", "VoltageSource"}, {"id", "narrow"}, {"params", {{"N", scalar}}}}.get<NarrowData>();
      success         *= rejects([&]
                         { EMT::parameter<uint8_t>(data, Parameter::N); },
                         "N");
    }
    SourceData wide;
    wide.parameters[Parameter::N]  = std::numeric_limits<size_t>::max();
    success                       *= EMT::parameter<size_t>(wide, Parameter::N) == std::numeric_limits<size_t>::max();
    wide.parameters[Parameter::N]  = static_cast<double>(std::numeric_limits<size_t>::max());
    success                       *= rejects([&]
                       { EMT::parameter<size_t>(wide, Parameter::N); });
    for (const auto value : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      wide.parameters[Parameter::omega]  = value;
      success                           *= rejects([&]
                         { EMT::parameter<double>(wide, Parameter::omega); },
                         "omega");
      wide.parameters[Parameter::E]      = EMT::ABCVector<double>{1, value, 3};
      success                           *= rejects([&]
                         { EMT::parameter<EMT::ABCVector<double>>(wide, Parameter::E); },
                         "E");
    }
    return success.report("Typed EMT numeric parameters");
  }

  Testing::TestOutcome strictModelData()
  {
    Testing::TestStatus success = true;
    const json          source  = {{"class", "VoltageSource"}, {"id", "source"}};
    for (const auto& patch : {
             json{{"params", {{"omeag", 3}}}}, json{{"params", {{"SIZE", 3}}}}, json{{"mon", {"unknown"}}}, json{{"mon", {"SIZE"}}}, json{{"mon", "ia"}}, json{{"inputs", {{"unknown", "bus.va"}}}}, json{{"inputs", {{"va", ""}}}}, json{{"outputs", {{"SIZE", "current"}}}}, json{{"outputs", {{"ia", ""}}}}, json{{"params", json::array({1, 2, 3})}}, json{{"unknown", 3}}, json{{"submodels", {{"unknown", json::object()}}}}})
    {
      auto input = source;
      input.update(patch);
      success *= rejects([&]
                         { input.get<SourceData>(); },
                         "source");
    }
    for (const auto& value : {json(nullptr), json("3"), json::array(), json::array({1, 2}), json::array({1, 2, 3, 4}), json::array({1, true, 3}), json::array({json::array({1, 2, 3}), json::array({1, 2}), json::array({1, 2, 3})}), json::array({1, json::array({2}), 3}), json(std::numeric_limits<double>::infinity())})
    {
      auto input            = source;
      input["params"]["E"]  = value;
      success              *= rejects([&]
                         { input.get<SourceData>(); },
                         "E");
    }
    const json model  = {{"header", {{"case_name", "strict"}, {"case_description", "input validation"}, {"case_comments", ""}}},
                         {"signals", {{{"id", "reference"}, {"value", 1}}}},
                         {"devices", json::array()}};
    success          *= model.get<EMT::SystemModelData<>>().signal[0].value == 1;
    for (const auto& patch : {
             json{{"unknown", true}}, json{{"devices", json::object()}}, json{{"signals", json::object()}}, json{{"signals", {{{"id", "reference"}, {"typo", 1}}}}}, json{{"monitors", {{{"format", "typo"}}}}}, json{{"monitors", {{{"format", "csv"}, {"typo", 1}}}}}, json{{"devices", {{{"class", "Container"}, {"id", "child"}, {"devices", json::array()}, {"typo", 1}}}}}})
    {
      auto input = model;
      input.update(patch);
      success *= rejects([&]
                         { input.get<EMT::SystemModelData<>>(); });
    }
    for (const auto& value : {json(true), json(std::numeric_limits<double>::infinity())})
    {
      const json fit  = {{"D", {{value, 0, 0}, {0, 1, 0}, {0, 0, 1}}}};
      success        *= rejects([&]
                         { fit.get<EMT::VectorFitData<double, size_t>>(); });
      success        *= rejects([&]
                         { fit.get<EMT::StateSpaceData<double, size_t>>(); });
    }
    const json unused  = {{"typo", 0}};
    success           *= rejects([&]
                       { unused.get<EMT::VectorFitData<double, size_t>>(); },
                       "typo");
    success           *= rejects([&]
                       { unused.get<EMT::StateSpaceData<double, size_t>>(); },
                       "typo");
    return success.report("Strict EMT model JSON fields and shapes");
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  results += numericParameters();
  results += strictModelData();
  return results.summary();
}
