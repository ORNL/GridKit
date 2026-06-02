/**
 * @file runForcedOscillationInjectionTests.cpp
 * @brief Unit tests for DynamicSimulation forced-oscillation study injections.
 */

#include <cmath>
#include <functional>
#include <iostream>
#include <optional>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"

namespace
{
  using namespace GridKit::PhasorDynamics;
  using AnalysisManager::Sundials::Ida;
  using GridKit::Testing::TestOutcome;
  using GridKit::Testing::TestStatus;

  using ModelData  = SystemModelData<>;
  using Trajectory = std::vector<std::vector<double>>;

  template <typename FuncT>
  bool throws(FuncT&& func)
  {
    try
    {
      func();
    }
    catch (const std::exception&)
    {
      return true;
    }
    return false;
  }

  void addSignal(ModelData& model, std::size_t id, const std::string& name)
  {
    ModelData::SignalDataT signal;
    signal.signal_id = id;
    signal.name      = name;
    model.signal.push_back(signal);
  }

  ModelData makeModel()
  {
    ModelData model;
    addSignal(model, 0, "speed");
    addSignal(model, 1, "pmech");
    addSignal(model, 2, "efd");

    ModelData::GenrouDataT gen;
    gen.device_class                                = "Genrou";
    gen.disambiguation_string                       = "GEN";
    gen.ports[ModelData::GenrouDataT::Ports::bus]   = 0;
    gen.ports[ModelData::GenrouDataT::Ports::speed] = 0;
    gen.ports[ModelData::GenrouDataT::Ports::pmech] = 1;
    gen.ports[ModelData::GenrouDataT::Ports::efd]   = 2;
    model.genrou.push_back(gen);

    ModelData::Tgov1DataT gov;
    gov.device_class                               = "Tgov1";
    gov.disambiguation_string                      = "GOV";
    gov.ports[ModelData::Tgov1DataT::Ports::speed] = 0;
    gov.ports[ModelData::Tgov1DataT::Ports::pmech] = 1;
    model.gov.push_back(gov);

    ModelData::Ieeet1DataT exciter;
    exciter.device_class                                = "Ieeet1";
    exciter.disambiguation_string                       = "EX";
    exciter.ports[ModelData::Ieeet1DataT::Ports::bus]   = 0;
    exciter.ports[ModelData::Ieeet1DataT::Ports::speed] = 0;
    exciter.ports[ModelData::Ieeet1DataT::Ports::efd]   = 2;
    model.exciter.push_back(exciter);

    return model;
  }

  ModelData makeTwoBusTgov1Model()
  {
    ModelData model;
    model.freq_base = 60.0;
    model.va_base   = 100.0e6;

    addSignal(model, 0, "speed");
    addSignal(model, 1, "pmech");

    auto& bus0    = model.bus.emplace_back();
    bus0.bus_id   = 0;
    bus0.bus_type = ModelData::BusDataT::BusType::DEFAULT;
    bus0.name     = "Bus 1";
    bus0.Vr0      = 0.9949877346411762;
    bus0.Vi0      = 0.09999703952427966;
    bus0.v_base   = 115.0e3;

    auto& bus1    = model.bus.emplace_back();
    bus1.bus_id   = 1;
    bus1.bus_type = ModelData::BusDataT::BusType::SLACK;
    bus1.name     = "Bus 2";
    bus1.Vr0      = 1.0;
    bus1.Vi0      = 0.0;
    bus1.v_base   = 115.0e3;

    using BranchParams                 = ModelData::BranchDataT::Parameters;
    using BranchPorts                  = ModelData::BranchDataT::Ports;
    auto& branch                       = model.branch.emplace_back();
    branch.device_class                = "Branch";
    branch.disambiguation_string       = "BR1";
    branch.ports[BranchPorts::bus1]    = 0;
    branch.ports[BranchPorts::bus2]    = 1;
    branch.parameters[BranchParams::R] = 0.0;
    branch.parameters[BranchParams::X] = 0.1;
    branch.parameters[BranchParams::G] = 0.0;
    branch.parameters[BranchParams::B] = 0.0;

    using GenParams                  = ModelData::GenrouDataT::Parameters;
    using GenPorts                   = ModelData::GenrouDataT::Ports;
    auto& gen                        = model.genrou.emplace_back();
    gen.device_class                 = "Genrou";
    gen.disambiguation_string        = "DV1";
    gen.ports[GenPorts::bus]         = 0;
    gen.ports[GenPorts::speed]       = 0;
    gen.ports[GenPorts::pmech]       = 1;
    gen.parameters[GenParams::p0]    = 1.0;
    gen.parameters[GenParams::q0]    = 0.05013;
    gen.parameters[GenParams::H]     = 3.0;
    gen.parameters[GenParams::D]     = 0.0;
    gen.parameters[GenParams::Ra]    = 0.0;
    gen.parameters[GenParams::Tdop]  = 7.0;
    gen.parameters[GenParams::Tdopp] = 0.04;
    gen.parameters[GenParams::Tqopp] = 0.05;
    gen.parameters[GenParams::Tqop]  = 0.75;
    gen.parameters[GenParams::Xd]    = 2.1;
    gen.parameters[GenParams::Xdp]   = 0.2;
    gen.parameters[GenParams::Xdpp]  = 0.18;
    gen.parameters[GenParams::Xq]    = 0.5;
    gen.parameters[GenParams::Xqp]   = 0.5;
    gen.parameters[GenParams::Xqpp]  = 0.18;
    gen.parameters[GenParams::Xl]    = 0.15;
    gen.parameters[GenParams::S10]   = 0.0;
    gen.parameters[GenParams::S12]   = 0.0;

    using GovParams                  = ModelData::Tgov1DataT::Parameters;
    using GovPorts                   = ModelData::Tgov1DataT::Ports;
    auto& gov                        = model.gov.emplace_back();
    gov.device_class                 = "Tgov1";
    gov.disambiguation_string        = "DV2";
    gov.ports[GovPorts::bus]         = 0;
    gov.ports[GovPorts::speed]       = 0;
    gov.ports[GovPorts::pmech]       = 1;
    gov.parameters[GovParams::R]     = 0.05;
    gov.parameters[GovParams::T1]    = 0.5;
    gov.parameters[GovParams::T2]    = 2.5;
    gov.parameters[GovParams::T3]    = 7.5;
    gov.parameters[GovParams::Pvmin] = 0.0;
    gov.parameters[GovParams::Pvmax] = 1.0;
    gov.parameters[GovParams::Dt]    = 0.0;

    return model;
  }

  ForcedOscillationInjection makeInjection()
  {
    ForcedOscillationInjection injection;
    injection.id     = "fo";
    injection.target = "Genrou:GEN.pmech";
    injection.mode   = "add";
    injection.params = {{"A", 0.02}, {"f", 0.7}};
    return injection;
  }

  void appendValues(std::vector<double>& sample, const std::vector<double>& values)
  {
    sample.insert(sample.end(), values.begin(), values.end());
  }

  std::vector<double> commonDynamicState(SystemModel<double, size_t>& system, bool has_fo)
  {
    std::vector<double> sample;
    appendValues(sample, system.getBus(0)->y());
    appendValues(sample, system.getBus(1)->y());

    const std::size_t gen_id = has_fo ? 2 : 1;
    const std::size_t gov_id = has_fo ? 3 : 2;

    auto* gen = dynamic_cast<Genrou<double, size_t>*>(system.getComponent(gen_id));
    auto* gov = dynamic_cast<Governor::Tgov1<double, size_t>*>(system.getComponent(gov_id));
    if (gen == nullptr || gov == nullptr)
    {
      return {};
    }

    appendValues(sample, gen->y());
    appendValues(sample, gov->y());
    return sample;
  }

  Trajectory runTrajectory(ModelData model, bool has_fo, double tmax, int steps)
  {
    SystemModel<double, size_t> system(model);
    system.allocate();

    Ida<double, size_t> ida(&system);
    ida.configureSimulation();
    ida.initializeSimulation(0.0, false);

    Trajectory trajectory;
    trajectory.push_back(commonDynamicState(system, has_fo));

    std::function<void(double)> callback = [&](double)
    {
      trajectory.push_back(commonDynamicState(system, has_fo));
    };
    ida.runSimulation(tmax, steps, std::optional<std::function<void(double)>>{callback});

    return trajectory;
  }

  bool trajectoriesMatch(const Trajectory& lhs, const Trajectory& rhs, double tolerance)
  {
    if (lhs.size() != rhs.size())
    {
      std::cout << "Trajectory length mismatch: " << lhs.size() << " vs " << rhs.size() << "\n";
      return false;
    }

    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
      if (lhs[i].size() != rhs[i].size() || lhs[i].empty())
      {
        std::cout << "State size mismatch at sample " << i << ": "
                  << lhs[i].size() << " vs " << rhs[i].size() << "\n";
        return false;
      }

      for (std::size_t j = 0; j < lhs[i].size(); ++j)
      {
        const double diff = std::abs(lhs[i][j] - rhs[i][j]);
        if (diff > tolerance)
        {
          std::cout << "Trajectory mismatch at sample " << i
                    << ", state " << j
                    << ": " << lhs[i][j]
                    << " vs " << rhs[i][j]
                    << " (diff=" << diff << ")\n";
          return false;
        }
      }
    }

    return true;
  }

  TestOutcome parser_shape()
  {
    TestStatus success = true;

    auto valid  = nlohmann::json::parse(R"json(
{
  "system_model_file": "case.json",
  "dt": 0.01,
  "tmax": 1.0,
  "forced_oscillations": []
}
)json");
    auto study  = StudyData(valid);
    success    *= study.isForcedOscillationStudy();
    success    *= study.events.empty();

    auto with_events  = nlohmann::json::parse(R"json(
{
  "system_model_file": "case.json",
  "dt": 0.01,
  "tmax": 1.0,
  "events": [],
  "forced_oscillations": []
}
)json");
    success          *= throws([&]()
                      { auto unused = StudyData(with_events); });

    auto with_unknown  = nlohmann::json::parse(R"json(
{
  "system_model_file": "case.json",
  "dt": 0.01,
  "tmax": 1.0,
  "forced_oscillations": [],
  "unknown": []
}
)json");
    success           *= throws([&]()
                      { auto unused = StudyData(with_unknown); });

    return success.report(__func__);
  }

  TestOutcome add_injection()
  {
    TestStatus success = true;

    auto base      = makeModel();
    auto injection = makeInjection();
    auto model     = applyInjections(base, {injection});

    using GenPorts = ModelData::GenrouDataT::Ports;
    using FoPorts  = ModelData::ForcedOscillationDataT::Ports;
    using FoMon    = ModelData::ForcedOscillationDataT::MonitorableVariables;

    success *= (base.signal.size() == 3);
    success *= (base.genrou[0].ports.at(GenPorts::pmech) == 1);
    success *= (model.signal.size() == 4);
    success *= (model.genrou[0].ports.at(GenPorts::pmech) == 3);
    success *= (model.forced_oscillation.size() == 1);
    success *= (model.forced_oscillation[0].ports.at(FoPorts::input) == 1);
    success *= (model.forced_oscillation[0].ports.at(FoPorts::output) == 3);
    success *= model.forced_oscillation[0].monitored_variables.contains(FoMon::in);
    success *= model.forced_oscillation[0].monitored_variables.contains(FoMon::force);
    success *= model.forced_oscillation[0].monitored_variables.contains(FoMon::out);
    success *= model.forced_oscillation[0].monitored_variables.contains(FoMon::active);

    return success.report(__func__);
  }

  TestOutcome drive_injection()
  {
    TestStatus success = true;

    auto injection   = makeInjection();
    injection.id     = "fo_vs";
    injection.target = "Ieeet1:EX.vs";
    injection.mode   = "drive";
    auto model       = applyInjections(makeModel(), {injection});

    using ExciterPorts = ModelData::Ieeet1DataT::Ports;
    using FoPorts      = ModelData::ForcedOscillationDataT::Ports;

    success *= (model.signal.size() == 4);
    success *= model.exciter[0].ports.contains(ExciterPorts::vs);
    success *= (model.exciter[0].ports.at(ExciterPorts::vs) == 3);
    success *= (model.forced_oscillation.size() == 1);
    success *= (!model.forced_oscillation[0].ports.contains(FoPorts::input));
    success *= (model.forced_oscillation[0].ports.at(FoPorts::output) == 3);

    return success.report(__func__);
  }

  TestOutcome signal_injection()
  {
    TestStatus success = true;

    auto injection   = makeInjection();
    injection.id     = "fo_signal";
    injection.target = "signal:9001";
    injection.mode   = "signal";
    auto model       = applyInjections(makeModel(), {injection});

    using FoPorts = ModelData::ForcedOscillationDataT::Ports;

    success *= (model.signal.size() == 4);
    success *= (model.signal.back().signal_id == 9001);
    success *= (model.forced_oscillation.size() == 1);
    success *= (model.forced_oscillation[0].ports.at(FoPorts::output) == 9001);

    return success.report(__func__);
  }

  TestOutcome invalid_injections()
  {
    TestStatus success = true;

    auto producer    = makeInjection();
    producer.target  = "Tgov1:GOV.pmech";
    success         *= throws([&]()
                      { auto unused = applyInjections(makeModel(), {producer}); });

    auto missing_port    = makeInjection();
    missing_port.target  = "Genrou:GEN.nope";
    success             *= throws([&]()
                      { auto unused = applyInjections(makeModel(), {missing_port}); });

    auto duplicate  = makeInjection();
    success        *= throws([&]()
                      { auto unused = applyInjections(makeModel(), {duplicate, duplicate}); });

    auto unknown_mode  = makeInjection();
    unknown_mode.mode  = "replace";
    success           *= throws([&]()
                      { auto unused = applyInjections(makeModel(), {unknown_mode}); });

    auto double_producer    = makeInjection();
    double_producer.target  = "signal:1";
    double_producer.mode    = "signal";
    success                *= throws([&]()
                      { auto unused = applyInjections(makeModel(), {double_producer}); });

    return success.report(__func__);
  }

  TestOutcome port_role_coverage()
  {
    TestStatus success  = true;
    success            *= forcedOscillationPortRoleTableComplete();
    return success.report(__func__);
  }

  TestOutcome max_step()
  {
    TestStatus success = true;

    auto dc    = makeInjection();
    dc.params  = {{"Bias", 1.0}};
    success   *= GridKit::Testing::isEqual(defaultForcedOscillationMaxStep({dc}, 0.01, 10.0), 0.01, 1.0e-12);

    auto chirp    = makeInjection();
    chirp.params  = {{"f", 1.0}, {"Kf", 1.0}, {"Ton", 1.0}};
    success      *= GridKit::Testing::isEqual(forcedOscillationMaxFrequency({chirp}, 3.0), 3.0, 1.0e-12);
    success      *= GridKit::Testing::isEqual(defaultForcedOscillationMaxStep({chirp}, 0.1, 3.0), 1.0 / 60.0, 1.0e-12);

    auto explicit_step  = nlohmann::json::parse(R"json(
{
  "system_model_file": "case.json",
  "dt": 0.01,
  "tmax": 1.0,
  "max_step": 0.005,
  "forced_oscillations": []
}
)json");
    auto study          = StudyData(explicit_step);
    success            *= study.max_step.has_value();
    success            *= GridKit::Testing::isEqual(*study.max_step, 0.005, 1.0e-12);

    return success.report(__func__);
  }

  TestOutcome delayed_add_identity_trajectory()
  {
    TestStatus success = true;

    auto injection   = makeInjection();
    injection.target = "Genrou:DV1.pmech";
    injection.params = {{"A", 0.02}, {"f", 0.7}, {"Ton", 10.0}, {"Tr", 1.0}};

    const auto   base_model = makeTwoBusTgov1Model();
    const auto   fo_model   = applyInjections(base_model, {injection});
    const double tmax       = 0.1;
    const int    steps      = 24;

    const auto base_trajectory = runTrajectory(base_model, false, tmax, steps);
    const auto fo_trajectory   = runTrajectory(fo_model, true, tmax, steps);

    success *= trajectoriesMatch(base_trajectory, fo_trajectory, 1.0e-6);

    return success.report(__func__);
  }

  TestOutcome zero_amplitude_identity_trajectory()
  {
    TestStatus success = true;

    auto injection   = makeInjection();
    injection.target = "Genrou:DV1.pmech";
    injection.params = {{"A", 0.0}, {"f", 0.7}, {"Ton", 0.0}};

    const auto   base_model = makeTwoBusTgov1Model();
    const auto   fo_model   = applyInjections(base_model, {injection});
    const double tmax       = 0.1;
    const int    steps      = 24;

    const auto base_trajectory = runTrajectory(base_model, false, tmax, steps);
    const auto fo_trajectory   = runTrajectory(fo_model, true, tmax, steps);

    success *= trajectoriesMatch(base_trajectory, fo_trajectory, 1.0e-6);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;

  result += parser_shape();
  result += add_injection();
  result += drive_injection();
  result += signal_injection();
  result += invalid_injections();
  result += port_role_coverage();
  result += max_step();
  result += delayed_add_identity_trajectory();
  result += zero_amplitude_identity_trajectory();

  return result.summary();
}
