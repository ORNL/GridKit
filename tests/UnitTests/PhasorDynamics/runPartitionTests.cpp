#include <exception>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "Partition.hpp"

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;
using DataT = SystemModelData<>;
using json  = nlohmann::json;

const json partition = json::parse(R"({
  "name": "root",
  "partitions": [
    {"name": "A", "internal": [], "external": [0, 1]},
    {"name": "group", "partitions": [
      {"name": "A", "internal": [2], "external": [0, 1]}
    ]}
  ]
})");

bool rejects(DataT data, const json& input)
{
  try
  {
    applyPartition(data, input);
  }
  catch (const std::exception&)
  {
    return true;
  }
  return false;
}

TestOutcome composition(DataT data)
{
  TestStatus success                                           = true;
  data.signal                                                  = {{"speed", 10}, {"pmech", 11}};
  data.genrou.at(1).signal_outputs[GenrouSignalOutputs::speed] = 10;
  data.genrou.at(1).signal_inputs[GenrouSignalInputs::pmech]   = 11;
  auto& governor                                               = data.gov.emplace_back();
  governor.signal_inputs[Governor::Tgov1SignalInputs::speed]   = 10;
  governor.signal_outputs[Governor::Tgov1SignalOutputs::pmech] = 11;
  data.monitor_sink.clear();
  auto composed                             = data;
  auto annotated                            = partition;
  annotated["description"]                  = "Nested runtime comparison";
  annotated["partitions"][0]["description"] = "Boundary-only partition";
  applyPartition(composed, annotated);

  const auto& a      = composed.subsystems.at(0);
  const auto& group  = composed.subsystems.at(1);
  const auto& b      = group.subsystems.at(0);
  success           *= composed.bus.size() == 2 && a.bus.empty() && b.bus.size() == 1;
  success           *= b.bus.at(0).bus_id == 2;
  success           *= a.branch.size() == 1 && b.branch.size() == 2;
  success           *= composed.branch.empty() && composed.loadz.empty() && b.loadz.size() == 1;
  success           *= a.level == 1 && b.level == 2 && a.id != b.id;
  success           *= b.freq_base == data.freq_base && b.va_base == data.va_base;
  success           *= composed.case_name == data.case_name;

  SystemModel<double, std::size_t> flat(data), nested(composed);
  success    *= flat.allocate() == 0 && nested.allocate() == 0;
  auto* leaf  = nested.getSubsystem(1)->getSubsystem(0);
  success    *= leaf->getBus(1) == nested.getBus(1);
  success    *= flat.size() == nested.size();
  success    *= nested.getBranch("BR_1_2") == leaf->getBranch("BR_1_2");
  success    *= leaf->getSignalNode(10)->linked() && leaf->getSignalNode(11)->linked();
  try
  {
    nested.getBranch("missing");
    success *= false;
  }
  catch (const std::out_of_range&)
  {
  }

  for (const bool fault : {false, true})
  {
    for (const bool in_service : {true, false, true})
    {
      for (auto* system : {&flat, &nested})
      {
        success *= system->initialize() == 0;
        success *= system->getBus(2)->setFault(fault, 0.0, 0.01) == 0;
        success *= system->getBranch("BR_1_2")->setInService(in_service) == 0;
        success *= system->evaluateResidual() == 0;
      }
      auto compare = [&](auto* original, auto* actual)
      {
        success *= original->size() == actual->size();
        for (std::size_t i = 0; i < original->size(); ++i)
        {
          success *= isEqual(original->y().getData()[i], actual->y().getData()[i]);
          success *= isEqual(flat.getResidual().getData()[original->getResidualIndex(i)],
                             nested.getResidual().getData()[actual->getResidualIndex(i)]);
        }
      };
      for (const auto& bus : data.bus)
      {
        compare(flat.getBus(bus.bus_id), nested.getBus(bus.bus_id));
      }
      compare(flat.getComponent(data.branch.size()), leaf->getComponent(b.branch.size()));
      for (std::size_t i = 0; i < data.genrou.size(); ++i)
      {
        const auto& subsystem_data = i == 0 ? a : b;
        auto*       subsystem      = i == 0 ? nested.getSubsystem(0) : leaf;
        compare(flat.getComponent(data.branch.size() + data.loadz.size() + i),
                subsystem->getComponent(subsystem_data.branch.size() + subsystem_data.loadz.size()));
      }
      compare(flat.getComponent(data.branch.size() + data.loadz.size() + data.genrou.size()),
              leaf->getComponent(b.branch.size() + b.loadz.size() + b.genrou.size()));
      success *= isEqual(flat.getSignalNode(11)->read(), leaf->getSignalNode(11)->read());
    }
  }
  return success.report(__func__);
}

TestOutcome signals(DataT data)
{
  TestStatus success                                           = true;
  data.signal                                                  = {{"pmech", 10}, {"pref", 11}, {"unused", 12}, {"source", 13}};
  data.genrou.at(0).signal_inputs[GenrouSignalInputs::pmech]   = 10;
  // Visit the unanchored controller first, so placement needs another pass.
  auto& upstream                                               = data.gov.emplace_back();
  upstream.signal_outputs[Governor::Tgov1SignalOutputs::pmech] = 11;
  upstream.signal_inputs[Governor::Tgov1SignalInputs::pref]    = 13;
  auto& governor                                               = data.gov.emplace_back();
  governor.signal_outputs[Governor::Tgov1SignalOutputs::pmech] = 10;
  governor.signal_inputs[Governor::Tgov1SignalInputs::pref]    = 11;
  auto& source                                                 = data.constant_source.emplace_back();
  source.signal_outputs[ConstantSignalSourceSignalOutputs::sr] = 13;

  auto composed = data;
  applyPartition(composed, partition);
  const auto& a  = composed.subsystems.at(0);
  success       *= a.genrou.size() == 1 && a.gov.size() == 2 && a.constant_source.size() == 1;
  success       *= a.signal.size() == 3 && composed.signal.size() == 1;
  success       *= composed.signal.at(0).signal_id == 12;

  auto bad                                                   = data;
  bad.genrou.at(1).signal_inputs[GenrouSignalInputs::pmech]  = 10;
  success                                                   *= rejects(bad, partition);
  bad                                                        = data;
  bad.signal.clear();
  success *= rejects(bad, partition);
  bad      = data;
  bad.constant_source.emplace_back();
  success                                            *= rejects(bad, partition);
  bad                                                 = data;
  bad.bus.at(1).signal_outputs[BusSignalOutputs::vr]  = 10;
  success                                            *= rejects(bad, partition);
  return success.report(__func__);
}

TestOutcome invalidPartitions(const DataT& data)
{
  TestStatus success = true;
  for (const auto& ids : {json::array({-1}), json::array({1.5}), json::array({99}), json::array({0, 0})})
  {
    auto bad                          = partition;
    bad["partitions"][0]["external"]  = ids;
    success                          *= rejects(data, bad);
  }
  auto bad                                           = partition;
  bad["partitions"][0]["internal"]                   = {2};
  success                                           *= rejects(data, bad);
  bad                                                = partition;
  bad["partitions"][0]["external"]                   = {0, 1, 2};
  success                                           *= rejects(data, bad);
  bad                                                = partition;
  bad["partitions"][1]["partitions"][0]["internal"]  = json::array();
  success                                           *= rejects(data, bad);
  bad                                                = partition;
  bad["partitions"][1]["partitions"][0]["external"]  = {0};
  success                                           *= rejects(data, bad);
  success                                           *= rejects(data, json{{"name", "empty"}, {"partitions", json::array()}});
  success                                           *= rejects(data, partition["partitions"][0]);
  for (const auto* key : {"internal", "external"})
  {
    bad       = partition;
    bad[key]  = json::array();
    success  *= rejects(data, bad);
  }

  auto composed = data;
  applyPartition(composed, partition);
  success *= rejects(composed, partition);
  return success.report(__func__);
}

int main(int argc, const char* argv[])
{
  if (argc != 2)
  {
    return 1;
  }
  const auto     data = parseSystemModelData(argv[1]);
  TestingResults result;
  result += composition(data);
  result += signals(data);
  result += invalidPartitions(data);
  return result.summary();
}
