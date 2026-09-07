#include <array>
#include <cmath>
#include <sstream>

#include <GridKit/Model/EMT/ComponentLibrary.hpp>
#include <GridKit/Model/EMT/DaeAnalysis.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using namespace GridKit::EMT;
  using System = SystemModel<double, size_t>;
  using Status = DaeAnalysis::Status;

  template <typename Action>
  bool rejects(Action action, const std::string& diagnostic)
  {
    try
    {
      action();
    }
    catch (const std::runtime_error& error)
    {
      return std::string(error.what()).find(diagnostic) != std::string::npos;
    }
    return false;
  }

  SystemModelData<double, size_t> circuit(double conductance, double capacitance)
  {
    std::istringstream stream(R"({
      "header": {"case_name":"DAE validation", "case_description":"", "case_comments":""},
      "devices": [{"class":"Container", "id":"plant", "devices":[
        {"class":"Bus", "id":"pcc", "outputs":{"va":"v_a"}},
        {"class":"VoltageSource", "id":"filter",
         "params":{"E":[1.0,1.0,1.0],"omega":1.0,
                   "Rs":[[2,0,0],[0,2,0],[0,0,2]],
                   "Ls":[[0.1,0,0],[0,0.1,0],[0,0,0.1]]},
         "inputs":{"va":"pcc.vb", "vb":"pcc.vc", "vc":"v_a"}}
      ], "signals":[{"id":"v_a"}]}]
    })");
    auto               data = parseSystemModelData(stream);
    if (conductance != 0.0 || capacitance != 0.0)
    {
      auto& Y = data.container[0].bus[0].shunts["shunt"];
      for (size_t p = 0; p < 3; ++p)
      {
        Y.D[p][p] = conductance;
        Y.E[p][p] = capacitance;
      }
    }
    return data;
  }

  bool matrixChecks()
  {
    // L i' + R i + v - e = 0; i - G v = 0.
    const JacobianEntries E{{0, 0, 0.1}};
    const JacobianEntries A{{0, 0, 2.0}, {0, 1, 1.0}, {1, 0, 1.0}};
    const auto            unsupported = analyzeDae(2, A, E);
    bool                  success     = unsupported.status == Status::structurally_singular
                   && unsupported.equations == std::vector<size_t>{1};
    auto resistive = A;
    resistive.push_back({1, 1, -0.2});
    success         &= analyzeDae(2, resistive, E).status == Status::regular;
    auto capacitive  = E;
    capacitive.push_back({1, 1, -0.01});
    success &= analyzeDae(2, A, capacitive).status == Status::regular;

    // Complete structural matching does not imply numerical independence.
    const JacobianEntries dependent{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
    success &= analyzeDae(2, {}, dependent).status == Status::numerically_singular;
    // Row 2 is the sum of rows 0 and 1. Rounded elimination leaves a tiny pivot.
    const JacobianEntries rounded{{0, 0, 1}, {0, 1, 2}, {0, 2, 3}, {1, 0, 4}, {1, 1, 5}, {1, 2, 6}, {2, 0, 5}, {2, 1, 7}, {2, 2, 9}};
    success &= analyzeDae(3, {}, rounded).status == Status::numerically_singular;
    // An incidental zero in F_y retains its structural slot and is a numerical failure.
    success &= analyzeDae(1, JacobianEntries{{0, 0, 0}}, {}).status == Status::numerically_singular;
    // Global accumulation can cancel a derivative coefficient.
    const JacobianEntries cancel{{0, 0, 1}, {0, 0, -1}};
    success &= analyzeDae(1, JacobianEntries{{0, 0, 2}}, cancel).differential.empty();
    success &= analyzeDae(0, {}, {}).status == Status::regular;
    return success;
  }

  bool systemChecks()
  {
    bool   success = true;
    System unsupported(circuit(0.0, 0.0));
    unsupported.allocate();
    unsupported.initialize();
    success &= rejects([&]
                       { unsupported.tagDifferentiable(); },
                       "Equation: plant.pcc.KCL[0]");

    for (const auto capacitance : {0.0, 0.01})
    {
      System system(circuit(0.2, capacitance));
      system.allocate();
      system.initialize();
      success      &= system.tagDifferentiable() == 0;
      auto& bus     = system.component("plant.pcc");
      auto& filter  = system.component("plant.filter");
      for (size_t p = 0; p < 3; ++p)
      {
        success &= bus.tag()[p] == (capacitance != 0.0);
        success &= system.tag()[bus.getVariableIndex(p)] == bus.tag()[p];
        success &= filter.tag()[3 + p];
      }

      const auto A = system.jacobianEntries(1.0, 0.0);
      const auto E = system.jacobianEntries(0.0, 1.0);
      using Matrix = std::map<std::pair<size_t, size_t>, double>;
      Matrix expected;
      for (const auto& e : A)
        expected[{e.row, e.column}] += e.value;
      for (const auto& e : E)
        expected[{e.row, e.column}] += 7.0 * e.value;
      for (const auto& e : system.jacobianEntries(1.0, 7.0))
        expected[{e.row, e.column}] -= e.value;
      for (const auto& [indices, value] : expected)
        success &= std::abs(value) < 1e-12;

      // Independent directional checks of each partial, including a permuted terminal.
      const size_t n = system.size();
      for (bool derivative : {false, true})
      {
        auto&               state = derivative ? system.yp() : system.y();
        std::vector<double> direction(n), actual(n, 0.0), before(n);
        for (size_t j = 0; j < n; ++j)
          direction[j] = 0.1 * static_cast<double>(j + 1);
        for (const auto& e : derivative ? E : A)
          actual[e.row] += e.value * direction[e.column];
        system.evaluateResidual();
        std::copy_n(system.getResidual().getData(), n, before.begin());
        constexpr double h = 1e-5;
        for (size_t j = 0; j < n; ++j)
          state.getData()[j] += h * direction[j];
        state.setDataUpdated();
        system.evaluateResidual();
        for (size_t j = 0; j < n; ++j)
          success &= std::abs((system.getResidual().getData()[j] - before[j]) / h - actual[j]) < 1e-8;
        for (size_t j = 0; j < n; ++j)
          state.getData()[j] -= h * direction[j];
        state.setDataUpdated();
      }
    }
    return success;
  }

  bool switchRestart()
  {
    auto  data                                        = circuit(0.0, 0.0);
    auto& plant                                       = data.container[0];
    auto& ground                                      = plant.bus.emplace_back();
    ground.id                                         = "ground";
    auto& source                                      = plant.voltage_source.emplace_back();
    source.id                                         = "grid";
    source.parameters[VoltageSourceParameters::omega] = 1.0;
    source.inputs                                     = {{VoltageSourceInputs::va, "ground.va"}, {VoltageSourceInputs::vb, "ground.vb"}, {VoltageSourceInputs::vc, "ground.vc"}};
    auto& sw                                          = plant.sw.emplace_back();
    sw.id                                             = "breaker";
    sw.inputs                                         = {{SwitchInputs::v1a, "pcc.va"}, {SwitchInputs::v1b, "pcc.vb"}, {SwitchInputs::v1c, "pcc.vc"}, {SwitchInputs::v2a, "ground.va"}, {SwitchInputs::v2b, "ground.vb"}, {SwitchInputs::v2c, "ground.vc"}};
    System system(data);
    system.allocate();
    system.initialize();
    system.getSwitch("plant.breaker")->setOpen(false);
    AnalysisManager::Sundials::Ida<double, size_t> ida(&system);
    ida.configureSimulation();
    ida.initializeSimulation(0.0);
    ida.runSimulation(0.001);
    system.getSwitch("plant.breaker")->setOpen(true);
    system.resetJacobianStructure();
    return rejects([&]
                   { ida.restartSimulation(0.001); },
                   "EMT DAE");
  }

  bool nonlinearSparsity()
  {
    std::istringstream stream(R"({
      "header":{"case_name":"Limiter", "case_description":"", "case_comments":""},
      "devices":[{"class":"Tgov1", "id":"governor",
                  "inputs":{"pref":"reference"}, "outputs":{"pmech":"power"}}],
      "signals":[{"id":"reference", "value":1.0}, {"id":"power"}]
    })");
    System             system(parseSystemModelData(stream));
    system.allocate();
    auto&                                  governor = system.component("governor");
    const auto                             valve    = governor.getVariableIndex(1);
    std::vector<std::pair<size_t, size_t>> previous;
    bool                                   success = true;
    for (const auto position : {10.0, 0.5, 10.0})
    {
      governor.y().getData()[1] = position;
      governor.y().setDataUpdated();
      const auto                             entries = system.jacobianEntries(1.0, 0.0);
      std::vector<std::pair<size_t, size_t>> pattern;
      double                                 slope = 0.0;
      for (const auto& entry : entries)
      {
        pattern.emplace_back(entry.row, entry.column);
        if (entry.row == valve && entry.column == valve)
          slope += entry.value;
      }
      success  &= slope == (position > 1.0 ? 0.0 : -2.0);
      success  &= previous.empty() || previous == pattern;
      previous  = std::move(pattern);
    }
    return success;
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  GridKit::Testing::TestStatus     status  = matrixChecks();
  results                                 += status.report("DAE matching, rank, and derivative classification");
  status                                   = systemChecks();
  results                                 += status.report("Assembled partials, nested aliases, and differential tags");
  status                                   = switchRestart();
  results                                 += status.report("Switching is validated before consistent initialization");
  status                                   = nonlinearSparsity();
  results                                 += status.report("Enzyme retains structural entries across limiter states");
  return results.summary();
}
