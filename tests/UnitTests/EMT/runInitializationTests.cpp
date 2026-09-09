#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit::EMT;
  using json   = nlohmann::json;
  using System = SystemModel<double, size_t>;
  using State  = std::map<std::string, std::map<std::string, double>>;

  auto model(const json& input)
  {
    std::istringstream stream(input.dump());
    return parseSystemModelData(stream);
  }

  // The 100 MVA, 13.8 kV, 50 MW equilibrium used by SystemTests.
  json machineCase()
  {
    return json::parse(R"({
      "header":{"case_name":"Initialization ownership", "case_description":"", "case_comments":""},
      "signals":[{"id":"speed"}, {"id":"pm"}, {"id":"efd"}],
      "devices":[
        {"class":"Bus", "id":"bus"},
        {"class":"Machine", "id":"machine",
         "params":{"N":3, "S":100e6, "V":13800.0, "f":60.0, "H":3.7, "F":0.0,
                   "Rs":0.003, "Ll":0.15, "Lmd":1.66, "Lmq":1.61, "L0":0.15,
                   "Rfd":0.0006, "Llfd":0.165, "R1d":0.0284, "Ll1d":0.1713,
                   "R1q":0.0062, "Ll1q":0.7252, "R2q":0.0237, "Ll2q":0.125,
                   "S10":0.1, "S12":0.5},
         "inputs":{"bus":"bus", "pm":"pm", "efd":"efd"}, "outputs":{"speed":"speed"}},
        {"class":"LoadZ", "id":"load", "inputs":{"bus":"bus"},
         "params":{"R":[[3.8088,0,0],[0,3.8088,0],[0,0,3.8088]]}},
        {"class":"Tgov1", "id":"governor", "inputs":{"speed":"speed"}, "outputs":{"pmech":"pm"}},
        {"class":"SexsPti", "id":"exciter", "inputs":{"bus":"bus"}, "outputs":{"efd":"efd"},
         "params":{"V":13800.0, "Ta":0.1, "Tb":0.5, "Te":0.2, "K":10.0, "Efdmin":-5.0, "Efdmax":5.0}}
      ]
    })");
  }

  State machineState(const std::string& prefix = "")
  {
    return {{prefix + "bus", {{"va", 11267.65281680262}, {"vb", -5633.82640840131}, {"vc", -5633.82640840131}}},
            {prefix + "machine", {{"ia", 2958.3209453903114}, {"ib", -1479.1604726951557}, {"ic", -1479.1604726951557}}}};
  }

  json nestedCase(bool reverse)
  {
    auto root    = machineCase();
    auto devices = root["devices"];
    auto plant   = json{{"class", "Container"}, {"id", "plant"}, {"inputs", {{"pm", "control.pm"}, {"efd", "control.efd"}}}, {"outputs", {{"speed", "speed"}}}, {"signals", json::array({{{"id", "speed"}}})}, {"devices", json::array({devices[0], devices[1], devices[2]})}};
    devices[4]["inputs"].erase("bus");
    for (const std::string phase : {"va", "vb", "vc"})
    {
      plant["outputs"][phase]     = "bus." + phase;
      devices[4]["inputs"][phase] = phase;
    }
    auto control = json{{"class", "Container"}, {"id", "control"}, {"inputs", {{"speed", "plant.speed"}, {"va", "plant.va"}, {"vb", "plant.vb"}, {"vc", "plant.vc"}}}, {"outputs", {{"pm", "pm"}, {"efd", "efd"}}}, {"signals", json::array({{{"id", "pm"}}, {{"id", "efd"}}})}, {"devices", json::array({devices[3], devices[4]})}};
    root.erase("signals");
    root["devices"] = json::array({plant, control});
    if (reverse)
    {
      std::reverse(root["devices"].begin(), root["devices"].end());
      for (auto& child : root["devices"])
        std::reverse(child["devices"].begin(), child["devices"].end());
    }
    return root;
  }

  bool near(double actual, double expected, double tolerance = 1e-10)
  {
    return std::isfinite(actual) && std::abs(actual - expected) < tolerance;
  }

  bool rejectsWithoutMutation(System& system, const State& state, const std::string& diagnostic, double omega = 0.0)
  {
    for (size_t n = 0; n < system.size(); ++n)
    {
      system.y().getData()[n]  = 1000.0 + static_cast<double>(n);
      system.yp().getData()[n] = -500.0 - static_cast<double>(n);
    }
    const std::vector<double> y(system.y().getData(), system.y().getData() + system.size());
    const std::vector<double> yp(system.yp().getData(), system.yp().getData() + system.size());
    bool                      rejected = false;
    try
    {
      system.initialize(state, omega);
    }
    catch (const std::invalid_argument& error)
    {
      rejected = std::string(error.what()).find(diagnostic) != std::string::npos;
    }
    return rejected && std::equal(y.begin(), y.end(), system.y().getData())
           && std::equal(yp.begin(), yp.end(), system.yp().getData());
  }

  bool nestedEquilibrium()
  {
    System flat(model(machineCase()));
    flat.allocate();
    bool success  = flat.initialize(machineState()) == 0;
    // Mechanical power includes the stator copper loss on the machine base.
    success      &= near(flat.signal("pm").read(), 0.5 + 0.003 * 0.5 * 0.5);
    success      &= near(flat.signal("efd").read(), 1.66 * flat.component("machine").y().getData()[12]);
    for (bool reverse : {false, true})
    {
      System nested(model(nestedCase(reverse)));
      nested.allocate();
      success &= nested.initialize(machineState("plant.")) == 0;
      for (const std::string& name : {"bus", "machine", "load", "governor", "exciter"})
      {
        const std::string prefix    = name == "governor" || name == "exciter" ? "control." : "plant.";
        auto&             expected  = flat.component(name);
        auto&             actual    = nested.component(prefix + name);
        success                    &= actual.size() == expected.size();
        for (size_t n = 0; n < expected.size(); ++n)
        {
          success &= near(actual.y().getData()[n], expected.y().getData()[n]);
          success &= near(actual.yp().getData()[n], expected.yp().getData()[n]);
        }
      }
      success     &= nested.evaluateResidual() == 0;
      double norm  = 0;
      for (size_t n = 0; n < nested.size(); ++n)
        norm += std::pow(nested.getResidual().getData()[n], 2);
      success &= std::sqrt(norm) < 1e-6;
    }
    return success;
  }

  bool prescribedConflicts()
  {
    bool success = true;
    for (const std::string& controller : {"governor", "exciter"})
    {
      System system(model(machineCase()));
      system.allocate();
      auto state                                                     = machineState();
      state[controller][controller == "governor" ? "pmech" : "efd"]  = 0.1;
      success                                                       &= rejectsWithoutMutation(system, state, controller);
    }
    auto input   = machineCase();
    auto second  = input["devices"][1];
    second["id"] = "second";
    second["inputs"].erase("efd");
    second.erase("outputs");
    input["devices"].push_back(second);
    System system(model(input));
    system.allocate();
    auto state      = machineState();
    state["second"] = state["machine"];
    for (auto& [name, value] : state["second"])
      value *= 0.5;
    success         &= rejectsWithoutMutation(system, state, "governor.pmech");
    // A shared producer remains valid when both consumers request the same value.
    state["second"]  = state["machine"];
    success         &= system.initialize(state) == 0;
    success         &= near(system.signal("pm").read(), 0.50075);
    return success;
  }

  bool constantOperatingPoint()
  {
    System reference(model(machineCase()));
    reference.allocate();
    reference.initialize(machineState());
    const double pm    = reference.signal("pm").read();
    const double efd   = reference.signal("efd").read();
    auto         input = machineCase();
    input["devices"].erase(4);
    input["devices"].erase(3);
    input["signals"][1]["value"] = pm;
    input["signals"][2]["value"] = efd;
    System matching(model(input));
    matching.allocate();
    bool success  = matching.initialize(machineState()) == 0;
    success      &= matching.signal("pm").constant() && matching.signal("pm").read() == pm;
    success      &= matching.signal("efd").constant() && matching.signal("efd").read() == efd;
    for (size_t signal : {size_t{1}, size_t{2}})
    {
      auto conflicting                        = input;
      conflicting["signals"][signal]["value"] = input["signals"][signal]["value"].get<double>() + 0.1;
      System system(model(conflicting));
      system.allocate();
      success &= rejectsWithoutMutation(system, machineState(), "constant " + conflicting["signals"][signal]["id"].get<std::string>());
    }
    return success;
  }

  bool invalidStates()
  {
    bool success = true;
    for (const auto& [path, output, value] :
         {std::tuple{"missing", "ia", 1.0}, std::tuple{"machine", "unknown", 1.0}, std::tuple{"bus", "va", std::numeric_limits<double>::infinity()}, std::tuple{"machine", "ia", std::numeric_limits<double>::quiet_NaN()}})
    {
      System system(model(machineCase()));
      system.allocate();
      auto state           = machineState();
      state[path][output]  = value;
      success             &= rejectsWithoutMutation(system, state, path);
    }
    auto input = machineCase();
    input["devices"].push_back({{"class", "Bus"}, {"id", "other"}});
    input["devices"].push_back({{"class", "Switch"}, {"id", "breaker"}, {"inputs", {{"bus1", "bus"}, {"bus2", "other"}}}});
    System system(model(input));
    system.allocate();
    auto state                = machineState();
    state["breaker"]["open"]  = 0.5;
    success                  &= rejectsWithoutMutation(system, state, "breaker");
    return success;
  }

  bool dependencyCycle()
  {
    System system(model(json::parse(R"({
      "header":{"case_name":"Initial dependency cycle", "case_description":"", "case_comments":""},
      "signals":[{"id":"a"}, {"id":"b"}],
      "devices":[
        {"class":"Tgov1", "id":"first", "inputs":{"speed":"b"}, "outputs":{"pmech":"a"}},
        {"class":"Tgov1", "id":"second", "inputs":{"speed":"a"}, "outputs":{"pmech":"b"}}
      ]
    })")));
    system.allocate();
    return rejectsWithoutMutation(system, {}, "Cyclic initialization dependency");
  }

  bool exciterReference()
  {
    auto input                             = machineCase();
    input["devices"][4]["class"]           = "Ieeet1";
    input["devices"][4]["params"]          = {{"V", 13800.0}, {"Ka", 10.0}, {"Ke", 1.0}, {"Vrmin", -5.0}, {"Vrmax", 5.0}, {"Se1", 0.0}, {"Se2", 0.0}};
    input["devices"][4]["inputs"]["speed"] = "speed";
    input["devices"][4]["inputs"]["vs"]    = "vs";
    input["signals"].push_back({{"id", "vs"}});
    input["devices"].push_back({{"class", "Ieeest"}, {"id", "stabilizer"}, {"params", {{"T6", 0.1}}}, {"inputs", {{"speed", "speed"}}}, {"outputs", {{"output", "vs"}}}});
    System reference(model(input));
    reference.allocate();
    bool         success = reference.initialize(machineState()) == 0;
    const double vref    = 1.0 + reference.signal("efd").read() / 10.0 - reference.signal("vs").read();
    input["signals"].push_back({{"id", "vref"}, {"value", vref}});
    input["devices"][4]["inputs"]["vref"] = "vref";
    System matching(model(input));
    matching.allocate();
    success &= matching.initialize(machineState()) == 0;
    matching.component("exciter").evaluateResidual();
    for (size_t n = 0; n < matching.component("exciter").size(); ++n)
      success &= near(matching.component("exciter").getResidual().getData()[n], 0.0, 1e-8);
    input["signals"].back()["value"] = vref + 0.1;
    System conflicting(model(input));
    conflicting.allocate();
    success &= rejectsWithoutMutation(conflicting, machineState(), "constant vref");
    return success;
  }

  json read(const std::filesystem::path& path)
  {
    std::ifstream stream(path);
    return json::parse(stream);
  }

  bool inverterEquilibrium()
  {
    const auto   directory = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().parent_path() / "cases/EMT/CurrentControl";
    const double omega     = 120 * std::acos(-1.0);
    bool         success   = true;
    for (const std::string name : {"GFL", "GFM"})
    {
      auto       input = read(directory / (name + ".case.json"));
      const auto data  = read(directory / (name + ".state.json"));
      State      state;
      for (const std::string section : {"buses", "devices"})
        for (const auto& [path, outputs] : data.at(section).items())
          state[path] = outputs.get<std::map<std::string, double>>();
      for (bool reverse : {false, true})
      {
        if (reverse)
          std::reverse(input["devices"].begin(), input["devices"].end());
        System system(model(input));
        system.allocate();
        success           &= system.initialize(state, omega) == 0;
        auto&      filter  = system.component("filter");
        auto&      inner   = system.component("current_control");
        auto&      outer   = system.component(name == "GFL" ? "power_control" : "voltage_control");
        // Independent dq circuit identity: the PI supplies the inverter-side copper drop.
        const auto params  = std::find_if(input["devices"].begin(), input["devices"].end(), [](const auto& device)
                                         { return device.at("id") == "filter"; })
                                ->at("params");
        const double resistance  = params.at("Rs")[0][0];
        const double id          = system.signal("id").read();
        const double iq          = system.signal("iq").read();
        success                 &= near(inner.y().getData()[0], resistance * id, 1e-8);
        success                 &= near(inner.y().getData()[1], resistance * iq, 1e-8);
        success                 &= near(inner.y().getData()[2], id, 1e-8);
        success                 &= near(inner.y().getData()[3], iq, 1e-8);
        success                 &= near(outer.y().getData()[2], id, 1e-8);
        success                 &= near(outer.y().getData()[3], iq, 1e-8);
        inner.evaluateResidual();
        outer.evaluateResidual();
        for (size_t n = 0; n < 2; ++n)
        {
          success &= near(inner.getResidual().getData()[n], 0.0, 1e-7);
          success &= near(outer.getResidual().getData()[n], 0.0, 1e-7);
        }
        // The physical LCL capacitor and grid inductor start on their sinusoidal orbit.
        filter.evaluateResidual();
        for (size_t n = 3; n < filter.size(); ++n)
          success &= near(filter.getResidual().getData()[n], 0.0, 1e-8);
        success                                                                   &= std::abs(filter.yp().getData()[4]) > 1.0;
        auto conflicting                                                           = state;
        conflicting[name == "GFL" ? "power_control" : "voltage_control"]["icmdd"]  = id + 1.0;
        success                                                                   &= rejectsWithoutMutation(system, conflicting, "icmdd", omega);
        conflicting                                                                = state;
        conflicting["filter"]["voa"]                                               = 1.0;
        success                                                                   &= rejectsWithoutMutation(system, conflicting, "filter.voa", omega);
      }
    }
    return success;
  }

  struct HistoryProbe : Container<double, size_t>
  {
    int                 resets = 0;
    double              bound  = std::numeric_limits<double>::infinity();
    std::vector<double> accepted;

    void resetHistory() override
    {
      ++resets;
      accepted.clear();
    }

    void acceptStep(double time) override
    {
      accepted.push_back(time);
    }

    double maximumStepSize() const override
    {
      return bound;
    }
  };

  bool nestedHistory()
  {
    Container<double, size_t> root;
    auto&                     nested   = root.add<Container<double, size_t>>("nested");
    auto&                     first    = nested.add<HistoryProbe>("first");
    auto&                     second   = root.add<HistoryProbe>("second");
    bool                      success  = std::isinf(root.maximumStepSize());
    first.bound                        = 0.1;
    second.bound                       = 0.025;
    success                           &= root.maximumStepSize() == 0.025;
    root.resetHistory();
    const std::vector<double> times{0.0, 0.05, 0.05, 0.1};
    for (double time : times)
      root.acceptStep(time);
    success &= first.resets == 1 && second.resets == 1;
    success &= first.accepted == times && second.accepted == times;
    root.resetHistory();
    success &= first.resets == 2 && second.resets == 2;
    success &= first.accepted.empty() && second.accepted.empty();
    for (double invalid : {0.0, -0.1, std::numeric_limits<double>::quiet_NaN()})
    {
      first.bound   = invalid;
      bool rejected = false;
      try
      {
        root.maximumStepSize();
      }
      catch (const std::invalid_argument&)
      {
        rejected = true;
      }
      success &= rejected;
    }
    return success;
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += GridKit::Testing::TestStatus(nestedEquilibrium()).report("nestedEquilibrium");
  result += GridKit::Testing::TestStatus(prescribedConflicts()).report("prescribedConflicts");
  result += GridKit::Testing::TestStatus(constantOperatingPoint()).report("constantOperatingPoint");
  result += GridKit::Testing::TestStatus(invalidStates()).report("invalidStates");
  result += GridKit::Testing::TestStatus(dependencyCycle()).report("dependencyCycle");
  result += GridKit::Testing::TestStatus(exciterReference()).report("exciterReference");
  result += GridKit::Testing::TestStatus(inverterEquilibrium()).report("inverterEquilibrium");
  result += GridKit::Testing::TestStatus(nestedHistory()).report("nestedHistory");
  return result.summary();
}
