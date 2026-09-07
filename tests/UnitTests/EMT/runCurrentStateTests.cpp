#include <array>
#include <limits>
#include <sstream>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

int main()
{
  using namespace GridKit;
  using json = nlohmann::json;
  Testing::TestingResults results;
  Testing::TestStatus     success = true;
  const auto              parse   = [](const json& value)
  {
    std::istringstream stream(value.dump());
    return Model::parseStateData(stream);
  };
  const auto rejects = [](auto&& action)
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

  const json         currents = {{"i12a", 4.0}, {"i12b", -3.0}, {"i12c", 0.0}, {"ia", -2.0}, {"ib", 1.5}, {"ic", 0.0}};
  const auto         state    = parse({{"devices", {{"network.device", currents}}}});
  std::ostringstream serialized;
  Model::writeStateData(state, serialized);
  success *= json::parse(serialized.str()).at("devices").at("network.device") == currents;
  std::istringstream round_trip(serialized.str());
  success *= Model::parseStateData(round_trip).devices.at("network.device").ib == 1.5;

  for (const auto* field : {"i12a", "i12b", "i12c", "ia", "ib", "ic"})
  {
    for (const auto& invalid : {json(true), json("1.0"), json::array({1.0}), json::object()})
    {
      success *= rejects([&]
                         { parse({{"devices", {{"device", {{field, invalid}}}}}}); });
    }
    const auto         empty = parse({{"devices", {{"device", {{field, nullptr}}}}}});
    std::ostringstream output;
    Model::writeStateData(empty, output);
    success *= json::parse(output.str()).at("devices").at("device").empty();
  }
  using CurrentMember = std::optional<double> Model::DeviceState::*;
  const std::array<CurrentMember, 6>          members{&Model::DeviceState::i12a, &Model::DeviceState::i12b, &Model::DeviceState::i12c, &Model::DeviceState::ia, &Model::DeviceState::ib, &Model::DeviceState::ic};
  for (const auto member : members)
  {
    for (const double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      Model::StateData invalid_state;
      invalid_state.devices["device"].*member  = invalid;
      success                                 *= rejects([&]
                         { std::ostringstream output; Model::writeStateData(invalid_state, output); });
    }
  }
  results += success.report("Current state parsing and serialization");

  const auto flat   = json::parse(R"({
    "header": {"case_name": "current state", "case_description": "state ingestion", "case_comments": ""},
    "devices": [
      {"class": "Bus", "id": "bus1"},
      {"class": "Bus", "id": "bus2"},
      {"class": "LineLumped", "id": "line", "params": {
        "dx": 1.0, "Rp": [[2,0,0],[0,2,0],[0,0,2]], "Lp": [[0.05,0,0],[0,0.05,0],[0,0,0.05]]},
        "inputs": {"bus1": "bus1", "bus2": "bus2"}},
      {"class": "LoadZ", "id": "load", "params": {
        "R": [[4,0,0],[0,4,0],[0,0,4]], "L": [[0.1,0,0],[0,0.1,0],[0,0,0.1]]},
        "inputs": {"bus": "bus1"}},
      {"class": "LoadZ", "id": "resistor", "params": {"R": [[5,0,0],[0,5,0],[0,0,5]]},
        "inputs": {"bus": "bus1"}}
    ]
  })");
  auto       nested = flat;
  nested["devices"].push_back({{"class", "Container"}, {"id", "network"}, {"devices", flat.at("devices")}});
  std::istringstream model_stream(nested.dump());
  auto               data = EMT::parseSystemModelData(model_stream);
  success                 = true;
  std::map<std::string, std::map<std::string, double>> initial;
  initial["bus1"] = initial["network.bus1"] = {{"va", 10.0}, {"vb", -5.0}, {"vc", -5.0}};
  initial["network.line"]                   = {{"i12a", 4.0}, {"i12b", -3.0}, {"i12c", 0.0}};
  initial["network.load"]                   = {{"ia", -2.0}, {"ib", 1.5}, {"ic", 0.0}};
  initial["network.resistor"]               = {{"ia", 99.0}};

  EMT::SystemModel<double, size_t> system(data);
  system.allocate();
  const std::array line_current{4.0, -3.0, 0.0};
  const std::array load_current{-2.0, 1.5, 0.0};
  for (int repeat = 0; repeat < 2; ++repeat)
  {
    system.component("network.line").y().setToConst(123.0);
    system.component("network.load").y().setToConst(123.0);
    success *= system.initialize(initial) == 0;
    for (size_t phase = 0; phase < 3; ++phase)
    {
      success *= system.component("network.line").y().getData()[phase] == line_current[phase];
      success *= system.component("network.load").y().getData()[phase] == load_current[phase];
      success *= system.component("line").y().getData()[phase] == 0.0;
      success *= system.component("load").y().getData()[phase] == 0.0;
    }
    success    *= system.component("network.resistor").y().getData()[0] == 99.0;
    success    *= system.component("network.resistor").y().getData()[1] == 1.0;
    success    *= system.component("resistor").y().getData()[0] == -2.0;
    auto& line  = system.component<EMT::LineLumped<double, size_t>>("network.line");
    success    *= line.size() == 3;
    for (size_t phase = 0; phase < 3; ++phase)
    {
      success *= line.outputSignal(static_cast<EMT::LineLumpedOutputs>(phase)).read() == line_current[phase];
      success *= line.outputSignal(static_cast<EMT::LineLumpedOutputs>(3 + phase)).read() == -line_current[phase];
    }
  }
  results += success.report("Nested current states, partial defaults, and repeated initialization");

  success = true;
  for (size_t phase = 0; phase < 3; ++phase)
  {
    for (const double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      auto invalid_state                                                = initial;
      invalid_state["network.line"][std::string("i12") + "abc"[phase]]  = invalid;
      success                                                          *= rejects([&]
                         { system.initialize(invalid_state); });
      invalid_state                                                     = initial;
      invalid_state["network.load"][std::string("i") + "abc"[phase]]    = invalid;
      success                                                          *= rejects([&]
                         { system.initialize(invalid_state); });
    }
  }
  results                                     += success.report("Nonfinite model current seeds are rejected");
  success                                      = true;
  auto explicit_case                           = flat;
  explicit_case["signals"]                     = json::array({{{"id", "v_a"}}, {{"id", "i_a"}}});
  explicit_case["devices"][0]["outputs"]       = {{"va", "v_a"}};
  explicit_case["devices"][1]["inputs"]        = {{"ia", "i_a"}};
  explicit_case["devices"][2]["inputs"]        = {{"v1a", "v_a"}, {"v1b", "bus1.vb"}, {"v1c", "bus1.vc"}, {"v2a", "bus2.va"}, {"v2b", "bus2.vb"}, {"v2c", "bus2.vc"}};
  explicit_case["devices"][2]["params"]["Cp"]  = {{0.001, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}};
  explicit_case["devices"][3]["inputs"]        = {{"va", "v_a"}, {"vb", "bus1.vb"}, {"vc", "bus1.vc"}};
  explicit_case["devices"][3]["outputs"]       = {{"ia", "i_a"}};
  std::istringstream explicit_stream(explicit_case.dump());
  auto               explicit_data  = EMT::parseSystemModelData(explicit_stream);
  success                          *= explicit_data.loadz[0].inputs.at(EMT::LoadZInputs::va) == "v_a";
  success                          *= data.loadz[0].inputs.at(EMT::LoadZInputs::vb) == "bus1.vb";
  EMT::SystemModel<double, size_t> explicit_system(explicit_data);
  success *= explicit_system.component<EMT::LineLumped<double, size_t>>("line")
                 .getSignals()
                 .getAttachedSignal<EMT::LineLumpedExternalVariables::V1A>()
             == &explicit_system.signal("v_a");
  success *= explicit_system.allocate() == 0;
  success *= explicit_system.initialize({{"bus1", {{"va", 10.0}, {"vb", -5.0}, {"vc", -5.0}}}, {"load", {{"ia", -2.0}}}}) == 0;
  explicit_system.tagDifferentiable();
  success *= explicit_system.component("bus1").tag()[0];
  success *= explicit_system.signal("v_a").read() == 10.0;
  success *= explicit_system.signal("i_a").read() == -2.0;
  // KCL includes both physical device contributions and an additional controlled injection.
  explicit_system.evaluateResidual();
  success                    *= explicit_system.component("bus1").getResidual().getData()[0] == -4.0;
  success                    *= explicit_system.component("bus2").getResidual().getData()[0] == -2.0;
  auto& bus2                  = explicit_system.component("bus2");
  success                    *= bus2.evaluateJacobian() == 0;
  auto* jacobian              = bus2.getCooJacobian();
  success                    *= jacobian != nullptr;
  double controlled_gradient  = 0.0;
  if (jacobian)
  {
    for (size_t k = 0; k < jacobian->getNnz(); ++k)
      if (jacobian->getRowData()[k] == bus2.getResidualIndex(0)
          && jacobian->getColData()[k] == explicit_system.component("load").getVariableIndex(0))
        controlled_gradient += jacobian->getValues()[k];
  }
  success *= controlled_gradient == 1.0;
  for (const auto* name : {"p", "SIZE"})
    success *= rejects([&]
                       { explicit_system.initialize({{"load", {{name, 1.0}}}}); });
  auto missing_kcl                   = explicit_case;
  missing_kcl["signals"][0]["value"] = 10.0;
  missing_kcl["devices"][0].erase("outputs");
  success                                 *= rejects([&]
                     {
                       std::istringstream stream(missing_kcl.dump());
                       EMT::SystemModel<double, size_t> invalid(EMT::parseSystemModelData(stream));
                       invalid.allocate(); invalid.initialize(); });
  auto duplicate                           = flat;
  duplicate["devices"][3]["inputs"]["va"]  = "bus1.va";
  success                                 *= rejects([&]
                     { std::istringstream stream(duplicate.dump()); EMT::parseSystemModelData(stream); });
  results                                 += success.report("Scalar phase inputs, bus shortcut expansion, output aliases, and KCL accumulation");

  success = true;
  EMT::Controller::Tgov1Data<double, size_t> governor_data;
  using GovernorParameter  = EMT::Controller::Tgov1Parameters;
  governor_data.parameters = {{GovernorParameter::R, 0.05}, {GovernorParameter::T1, 0.1}, {GovernorParameter::T2, 0.2}, {GovernorParameter::T3, 0.5}, {GovernorParameter::Pvmin, 0.0}, {GovernorParameter::Pvmax, 1.0}, {GovernorParameter::Dt, 0.0}};
  explicit_data.gov.push_back(governor_data);
  explicit_data.gov.back().id = "governor";
  EMT::Controller::Tgov1<double, size_t> governor(explicit_data.gov.back());
  EMT::Signal<double, size_t>            pmech;
  governor.getSignals().template assignSignal<EMT::Controller::Tgov1InternalVariables::PM>(&pmech);
  success *= governor.allocate() == 0 && governor.initialize({{EMT::Controller::Tgov1Outputs::pmech, 0.7}}) == 0;
  governor.evaluateResidual();
  for (size_t row = 0; row < 3; ++row)
    success *= std::abs(governor.getResidual().getData()[row]) < 1e-12;
  success *= governor.y().getData()[2] == 0.7;
  results += success.report("Controller output state back-calculates its operating point");

  success = true;
  EMT::Controller::PwmData<double, size_t> pwm_data;
  using PwmParameter  = EMT::Controller::PwmParameters;
  pwm_data.parameters = {{PwmParameter::M, 0.8}, {PwmParameter::fm, 60.0}, {PwmParameter::fc, 900.0}, {PwmParameter::alignment, 0.5}};
  EMT::Controller::Pwm<double, size_t> pwm(pwm_data);
  pwm.allocate();
  success *= pwm.initialize({{EMT::Controller::PwmOutputs::sa, static_cast<double>(pwm.output(0))}}) == 0;
  success *= rejects([&]
                     { pwm.initialize({{EMT::Controller::PwmOutputs::sa, 20.0}}); });
  results += success.report("Computed outputs reject incompatible initial values");
  success  = true;
  EMT::SystemModelData<double, size_t> fitted_data;
  auto&                                fitted_bus = fitted_data.bus.emplace_back();
  fitted_bus.id                                   = "bus";
  auto& source                                    = fitted_data.voltage_source.emplace_back();
  using SourceInput                               = EMT::VoltageSourceInputs;
  using SourceOutput                              = EMT::VoltageSourceOutputs;
  source.id                                       = "source";
  source.inputs                                   = {{SourceInput::va, "bus.va"}, {SourceInput::vb, "bus.vb"}, {SourceInput::vc, "bus.vc"}};
  source.parameters                               = {{EMT::VoltageSourceParameters::E, EMT::ABCVector<double>{20.0, 20.0, 20.0}},
                                                     {EMT::VoltageSourceParameters::omega, 20.0}};
  source.Y.emplace();
  source.Y->poles = {{-3.0, 0.0}};
  source.Y->residues.resize(1);
  for (size_t n = 0; n < 3; ++n)
  {
    source.Y->D[n][n]           = 2.0;
    source.Y->residues[0][n][n] = 4.0;
  }
  const double                     expected_current = 2.0 * (std::sqrt(2.0) * 20.0 - 10.0);
  EMT::SystemModel<double, size_t> fitted_system(fitted_data);
  success                                   *= fitted_system.allocate() == 0;
  success                                   *= fitted_system.initialize({{"bus", {{"va", 10.0}, {"vb", -5.0}, {"vc", -5.0}}}, {"source", {{"ia", expected_current}}}}) == 0;
  auto&                       fitted_source  = fitted_system.component<EMT::VoltageSource<double, size_t>>("source");
  EMT::Signal<double, size_t> fitted_current;
  fitted_source.assignOutput(SourceOutput::ia, &fitted_current);
  success *= std::abs(fitted_current.read() - expected_current) < 1e-12;
  EMT::Signal<double, size_t>::GradientT gradient;
  fitted_current.appendGradient(gradient);
  success *= gradient.size() == 2;
  for (const auto& [column, coefficient] : gradient)
  {
    auto&        value     = fitted_system.y().getData()[column];
    const double original  = value;
    value                  = original + 1e-3;
    const double plus      = fitted_current.read();
    value                  = original - 1e-3;
    const double minus     = fitted_current.read();
    value                  = original;
    success               *= std::abs((plus - minus) / 2e-3 - coefficient) < 1e-9;
  }
  success *= rejects([&]
                     { fitted_source.initialize({{SourceOutput::ia, 1.0}}); });
  results += success.report("Rational current outputs preserve voltage states and include memory gradients");
  return results.summary();
}
