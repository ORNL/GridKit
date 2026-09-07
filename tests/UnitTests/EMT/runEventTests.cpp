#include <cmath>
#include <limits>
#include <sstream>

#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <application/EMT/EventSchedule.hpp>

namespace
{
  using namespace GridKit::EMT;
  using System   = SystemModel<double, size_t>;
  using Ida      = AnalysisManager::Sundials::Ida<double, size_t>;
  using Schedule = EventSchedule<double, size_t>;

  bool near(double actual, double expected, double tolerance = 1e-7)
  {
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance;
  }

  bool rejects(auto&& action)
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
  }

  auto model(const json& input)
  {
    std::istringstream stream(input.dump());
    return parseSystemModelData(stream);
  }

  auto governorModel()
  {
    return model(json::parse(R"({
      "header":{"case_name":"Reference step", "case_description":"", "case_comments":""},
      "devices":[{"class":"Container", "id":"plant",
        "signals":[{"id":"reference", "value":0.2}, {"id":"power"}],
        "devices":[{"class":"Tgov1", "id":"governor",
          "params":{"R":1.0, "T1":0.5, "T2":0.25, "T3":0.25,
                    "Pvmin":-10.0, "Pvmax":10.0, "Dt":0.0},
          "inputs":{"pref":"reference"}, "outputs":{"pmech":"power"}}]}]
    })"));
  }

  StudyData study()
  {
    auto result  = json{{"system_model_file", "unused.json"}, {"tmax", 1.0}, {"dt_monitor", 0.05}}.get<StudyData>();
    result.state = {{"plant.governor", {{"pmech", 0.2}}}};
    return result;
  }

  bool constantBinding()
  {
    using SignalT = Signal<double, size_t>;
    SignalT signal("reference");
    signal.bindConstant(0.2);
    bool success = signal.constant() && signal.computed() && signal.linked()
                   && !signal.derivativeLinked() && !signal.residualLinked();
    SignalT::GradientT gradient;
    signal.appendGradient(gradient);
    success &= gradient.empty();
    signal.setConstantValue(0.8);
    success &= signal.read() == 0.8;
    success &= rejects([&]
                       { signal.init(0.1); });
    success &= rejects([&]
                       { signal.bindConstant(0.1); });
    success &= rejects([&]
                       { signal.claimProducer(); });
    for (double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      success &= rejects([&]
                         { signal.setConstantValue(invalid); });
      success &= signal.read() == 0.8;
      SignalT unbound;
      success &= rejects([&]
                         { unbound.bindConstant(invalid); });
      success &= !unbound.hasProducer();
    }
    double value = 0.5;
    size_t index = 0;
    signal.set(&value, &index);
    success &= !signal.constant() && signal.read() == value;
    success &= rejects([&]
                       { signal.setConstantValue(1.0); });
    success &= rejects([&]
                       { signal.bindConstant(1.0); });
    SignalT expression;
    expression.bindConstant(0.5);
    expression.setComputed([]
                           { return 2.0; },
                           [](auto&, double) {});
    success &= !expression.constant() && expression.read() == 2.0;
    success &= rejects([&]
                       { expression.setConstantValue(1.0); });
    SignalT unbound;
    success &= rejects([&]
                       { unbound.setConstantValue(1.0); });
    return success;
  }

  bool referenceStep()
  {
    auto data               = study();
    data.consistent_ic_type = AnalysisManager::Sundials::IdaConsistentICType::Y;
    data.events             = {{0.25, SignalStep{"plant.reference", 0.6}},
                               {0.25, SignalStep{"plant.reference", 0.8}}};
    System   system(governorModel());
    Schedule events(system, data);
    system.allocate();
    system.initialize(data.state);
    Ida ida(&system);
    ida.setTolerance(1e-10, 1e-12);
    ida.setConsistentICType(data.consistent_ic_type);
    events.configure(ida);
    auto& governor = system.component("plant.governor");

    struct Sample
    {
      double time, value, derivative, reference;
    };

    std::vector<Sample> samples;
    events.run(ida, [&](double time)
               { samples.push_back({time, governor.y().getData()[1], governor.yp().getData()[1], system.signal("plant.reference").read()}); });
    bool                success = true;
    std::vector<Sample> at_event;
    for (const auto& sample : samples)
    {
      // T2=T3 makes the turbine output follow the valve's first-order response.
      const double expected  = sample.time <= 0.25 ? 0.2 : 0.8 - 0.6 * std::exp(-(sample.time - 0.25) / 0.5);
      success               &= near(sample.value, expected);
      if (sample.time == 0.25)
        at_event.push_back(sample);
    }
    success &= at_event.size() == 2;
    if (at_event.size() == 2)
    {
      success &= near(at_event[0].value, at_event[1].value, 1e-11);
      success &= near(at_event[0].derivative, 0.0, 1e-10);
      success &= near(at_event[1].derivative, 1.2, 1e-8);
      success &= at_event[0].reference == 0.2 && at_event[1].reference == 0.8;
    }
    success &= near(system.signal("plant.power").read(), 0.8 - 0.6 * std::exp(-1.5));
    return success;
  }

  bool constantReferences()
  {
    System system(model(json::parse(R"({
      "header":{"case_name":"Constant references", "case_description":"", "case_comments":""},
      "signals":[{"id":"pref", "value":0.3}, {"id":"vref", "value":1.1},
                 {"id":"steam_power"}, {"id":"gas_power"}, {"id":"sexs_field"}, {"id":"ieee_field"}],
      "devices":[{"class":"Bus", "id":"bus"},
        {"class":"VoltageSource", "id":"source", "params":{"omega":1.0}, "inputs":{"bus":"bus"}},
        {"class":"Tgov1", "id":"steam", "inputs":{"pref":"pref"}, "outputs":{"pmech":"steam_power"}},
        {"class":"GastPti", "id":"gas", "params":{"S":100000000.0},
         "inputs":{"pref":"pref"}, "outputs":{"pmech":"gas_power"}},
        {"class":"SexsPti", "id":"sexs",
         "params":{"V":100.0, "Ta":0.1, "Tb":0.5, "Te":0.2, "K":10.0, "Efdmin":-5.0, "Efdmax":5.0},
         "inputs":{"bus":"bus", "vref":"vref"}, "outputs":{"efd":"sexs_field"}},
        {"class":"Ieeet1", "id":"ieee", "params":{"V":100.0, "Ke":1.0},
         "inputs":{"bus":"bus", "vref":"vref"}, "outputs":{"efd":"ieee_field"}}]
    })")));
    system.allocate();
    system.initialize({{"steam", {{"pmech", 0.2}}}, {"gas", {{"pmech", 0.2}}}, {"sexs", {{"efd", 0.2}}}, {"ieee", {{"efd", 0.2}}}});
    return system.signal("pref").read() == 0.3 && system.signal("vref").read() == 1.1
           && near(system.signal("steam_power").read(), 0.2)
           && near(system.signal("gas_power").read(), 0.2)
           && near(system.signal("sexs_field").read(), 0.2)
           && near(system.signal("ieee_field").read(), 0.2);
  }

  bool timeZero()
  {
    auto data   = study();
    data.tmax   = 0.0;
    data.events = {{0.0, SignalStep{"plant.reference", 0.4}}};
    System   system(governorModel());
    Schedule events(system, data);
    system.allocate();
    system.initialize(data.state);
    bool success = system.signal("plant.reference").read() == 0.2;
    Ida  ida(&system);
    ida.setTolerance(1e-10, 1e-12);
    events.configure(ida);
    success        &= system.signal("plant.reference").read() == 0.4;
    size_t samples  = 0;
    events.run(ida, [&](double time)
               {
                 ++samples;
                 auto& governor = system.component("plant.governor");
                 success &= time == 0.0 && near(governor.y().getData()[1], 0.2, 1e-11);
                 success &= near(governor.yp().getData()[1], 0.4, 1e-9); });
    return success && samples == 1;
  }

  bool timeZeroSwitch()
  {
    auto data   = study();
    data.state  = {{"breaker", {{"open", 1.0}}}};
    data.tmax   = 0.001;
    data.events = {{0.0, SwitchEvent{"breaker", false}}};
    System   system(model(json::parse(R"({
      "header":{"case_name":"Initial switch", "case_description":"", "case_comments":""},
      "devices":[{"class":"Bus", "id":"grid"}, {"class":"Bus", "id":"island"},
        {"class":"VoltageSource", "id":"source", "params":{"omega":1.0}, "inputs":{"bus":"grid"}},
        {"class":"Switch", "id":"breaker", "inputs":{"bus1":"grid", "bus2":"island"}}]
    })")));
    Schedule events(system, data);
    system.allocate();
    system.initialize(data.state);
    // The open breaker leaves a free bus voltage. The time-zero closure must
    // precede the first DAE check in configureSimulation().
    Ida ida(&system);
    events.configure(ida);
    events.run(ida);
    return true;
  }

  bool invalidSchedules()
  {
    System system(governorModel());
    bool   success = true;
    for (const auto& action : std::vector<std::variant<SwitchEvent, SignalStep>>{
             SignalStep{"missing", 0.3}, SignalStep{"plant.power", 0.3}, SignalStep{"plant.reference", std::numeric_limits<double>::infinity()}, SwitchEvent{"missing", true}, SwitchEvent{"plant.governor", true}})
    {
      auto data    = study();
      data.events  = {{0.0, SignalStep{"plant.reference", 0.4}}, {0.5, action}};
      success     &= rejects([&]
                         { Schedule events(system, data); });
      success     &= system.signal("plant.reference").read() == 0.2;
    }
    for (double time : {-0.1, 1.1, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      auto data    = study();
      data.events  = {{time, SignalStep{"plant.reference", 0.4}}};
      success     &= rejects([&]
                         { Schedule events(system, data); });
    }
    auto data    = study();
    data.events  = {{0.5, SignalStep{"plant.reference", 0.4}}, {0.25, SignalStep{"plant.reference", 0.6}}};
    success     &= rejects([&]
                       { Schedule events(system, data); });
    return success;
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  GridKit::Testing::TestStatus     status  = constantBinding();
  results                                 += status.report("Declared constant ownership and rebinding");
  status                                   = referenceStep();
  results                                 += status.report("Grouped reference steps preserve states and match a first-order response");
  status                                   = constantReferences();
  results                                 += status.report("Controller initialization preserves supplied constant references");
  status                                   = timeZero();
  results                                 += status.report("Time-zero reference step and zero-duration study");
  status                                   = timeZeroSwitch();
  results                                 += status.report("Time-zero switch precedes DAE validation");
  status                                   = invalidSchedules();
  results                                 += status.report("Invalid schedules fail before applying any event");
  return results.summary();
}
