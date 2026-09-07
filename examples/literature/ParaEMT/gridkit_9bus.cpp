// Case-specific initialization and disturbance driver; component equations are unchanged.
#include <chrono>
#include <iomanip>
#include <map>

#include <GridKit/Model/EMT/SystemModel.hpp>

#include "GeneratorTerminalConstraint.hpp"
#include <application/EMT/AnalysisUtilities.hpp>
#include <application/EMT/StateMonitor.hpp>

using json   = nlohmann::json;
using System = GridKit::EMT::SystemModel<double, size_t>;
using namespace AnalysisManager::Sundials;

// The application normally starts passive inductors at zero current. This
// experiment starts at the supplied balanced power flow. The phasor seeds
// include currents AND derivatives; neither a warm-up nor frozen controls is used.
class NineBus : public System
{
  json                       seed_;
  std::array<std::string, 3> transformers_{"xfmr_4_1", "xfmr_7_2", "xfmr_9_3"};
  std::array<double, 3>      current_base_{};
  std::array<double, 3>      pref_{};
  bool                       trip_case_, tripped_{false};
  std::array<double, 2>      open_d_{}, open_q_{};
  size_t                     constant_index_{GridKit::INVALID_INDEX<size_t>};

public:
  NineBus(const GridKit::EMT::StudyData& study, const std::filesystem::path& seed, bool trip)
    : System(study.model_data), seed_(json::parse(std::ifstream(seed))), trip_case_(trip)
  {
    if (trip)
      seed_["machine_bus_1"] = seed_.at("bus_1");
    // Writable setpoints receive each governor's operating-point seed.
    for (size_t i = 0; i < pref_.size(); ++i)
    {
      getSignal("pref_" + std::to_string(i + 1))->set(&pref_[i], &constant_index_);
      using P                  = GridKit::EMT::MachineParameters;
      const auto& machine_data = study.model_data.machine.at(i);
      if (machine_data.id != "gen_" + std::to_string(i + 1))
        throw std::invalid_argument("Unexpected machine ordering");
      const auto& parameters = machine_data.parameters;
      if (std::get<double>(parameters.at(P::S10)) != 0 || std::get<double>(parameters.at(P::S12)) != 0)
        throw std::invalid_argument("Terminal index reduction requires the unsaturated case");
      const auto transformer = std::find_if(study.model_data.line_lumped.begin(), study.model_data.line_lumped.end(), [&](const auto& line)
                                            { return line.id == transformers_[i]; });
      if (transformer == study.model_data.line_lumped.end())
        throw std::invalid_argument("Missing generator transformer");
      for (auto parameter : {GridKit::EMT::LineLumpedParameters::Cp, GridKit::EMT::LineLumpedParameters::Gp})
        for (const auto& row : std::get<GridKit::EMT::ABCMatrix<double>>(transformer->parameters.at(parameter)))
          for (double value : row)
            if (value != 0)
              throw std::invalid_argument("Terminal index reduction requires zero transformer shunts");
      current_base_[i]        = std::sqrt(2.0 / 3.0) * std::get<double>(parameters.at(P::S)) / std::get<double>(parameters.at(P::V));
      const auto coefficients = json::parse(std::ifstream(seed.parent_path() / "index_reduction.json")).at(std::to_string(i + 1));
      add<GeneratorTerminalConstraint>("terminal_constraint_" + std::to_string(i + 1),
                                       component("bus_" + std::to_string(i + 1)),
                                       component("gen_" + std::to_string(i + 1)),
                                       component(transformers_[i]),
                                       coefficients.at("d").get<std::array<double, 3>>(),
                                       coefficients.at("q").get<std::array<double, 3>>(),
                                       current_base_[i],
                                       coefficients.at("L0").get<double>(),
                                       coefficients.at("omega").get<double>(),
                                       trip && i == 0 ? &component("machine_bus_1") : nullptr);
      if (i == 0)
      {
        auto value = [&](P key)
        { return std::get<double>(parameters.at(key)); };
        auto open_flux = [](double lm, double l1, double l2)
        {
          const double determinant = (lm + l1) * (lm + l2) - lm * lm;
          return std::array<double, 2>{lm * l2 / determinant, lm * l1 / determinant};
        };
        open_d_ = open_flux(value(P::Lmd), value(P::Llfd), value(P::Ll1d));
        open_q_ = open_flux(value(P::Lmq), value(P::Ll1q), value(P::Ll2q));
      }
    }
  }

  int initialize() override
  {
    // getDefaultInitialCondition is IDA's public model-to-solver copy path.
    // After the explicit opening projection, its state is already initialized.
    if (tripped_)
      return 0;
    const int status = System::initialize();
    if (status != 0)
      return status;
    for (const auto& [name, data] : seed_.items())
    {
      auto& device = component(name);
      if (data.at("y").size() != device.size() || data.at("yp").size() != device.size())
        throw std::runtime_error("Initial state size mismatch for " + name);
      for (size_t i = 0; i < device.size(); ++i)
      {
        device.y().getData()[i]  = data.at("y").at(i).get<double>();
        device.yp().getData()[i] = data.at("yp").at(i).get<double>();
      }
    }
    for (size_t i = 0; i < 3; ++i)
    {
      auto& generator   = component("gen_" + std::to_string(i + 1));
      auto& transformer = component(transformers_[i]);
      for (size_t phase = 0; phase < 3; ++phase)
        generator.yp().getData()[21 + phase] = -transformer.yp().getData()[phase] / current_base_[i];
    }
    y().setDataUpdated();
    yp().setDataUpdated();
    return 0;
  }

  double generatorCurrentMismatch()
  {
    double maximum = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      auto& generator   = component("gen_" + std::to_string(i + 1));
      auto& transformer = component(transformers_[i]);
      for (size_t phase = 0; phase < 3; ++phase)
      {
        const double machine_current = current_base_[i] * generator.y().getData()[21 + phase];
        maximum                      = std::max(maximum, std::abs((tripped_ && i == 0 ? 0.0 : machine_current) + transformer.y().getData()[phase] + transformer.y().getData()[6 + phase]));
        if (tripped_ && i == 0)
          maximum = std::max(maximum, std::abs(machine_current));
      }
    }
    return maximum;
  }

  json tripGenerator()
  {
    if (!trip_case_ || tripped_)
      throw std::logic_error("Unexpected trip");
    auto& machine     = component("gen_1");
    auto& transformer = component("xfmr_4_1");
    auto* y           = machine.y().getData();
    json  changes     = json::array();
    auto  project     = [&](auto& device, size_t index, double value)
    {
      changes.push_back({{"index", device.getVariableIndex(index)}, {"before", device.y().getData()[index]}, {"after", value}});
      device.y().getData()[index] = value;
    };
    // Ideal current interruption: stator currents and transformer currents
    // become zero. Preserve rotor winding fluxes, rotor motion, controls and
    // all other network differential states. Terminal voltage has an impulse;
    // this is a pre/post event solution, not a finite breaker-arc waveform.
    project(machine, 2, open_d_[0] * y[5] + open_d_[1] * y[6]);
    project(machine, 3, open_q_[0] * y[7] + open_q_[1] * y[8]);
    project(machine, 4, 0.0);
    for (size_t phase = 0; phase < 3; ++phase)
      project(transformer, phase, 0.0);
    dynamic_cast<GeneratorTerminalConstraint&>(component("terminal_constraint_1")).open();
    tripped_ = true;
    this->y().setDataUpdated();
    yp().setDataUpdated();
    return changes;
  }
};

// Independently check the assembled residual/Jacobian, including every
// derivative column introduced by the terminal constraint, at two run states.
double checkJacobian(NineBus& sys, double time)
{
  const size_t n     = sys.size();
  const double alpha = 123.45;
  sys.updateTime(time, alpha);
  sys.evaluateResidual();
  sys.evaluateJacobian();
  auto*                                       jac = sys.getCsrJacobian();
  std::map<std::pair<size_t, size_t>, double> entries;
  for (size_t row = 0; row < n; ++row)
    for (size_t k = jac->getRowData()[row]; k < jac->getRowData()[row + 1]; ++k)
      entries[{row, jac->getColData()[k]}] = jac->getValues()[k];
  double maximum = 0;
  for (size_t column = 0; column < n; ++column)
  {
    double*      y     = sys.y().getData();
    double*      yp    = sys.yp().getData();
    const double value = y[column], derivative = yp[column];
    // Keep phase-angle perturbations small even after many electrical cycles.
    const double h = 1e-6;
    y[column]      = value + h;
    yp[column]     = derivative + alpha * h;
    sys.evaluateResidual();
    std::vector<double> plus(sys.getResidual().getData(), sys.getResidual().getData() + n);
    y[column]  = value - h;
    yp[column] = derivative - alpha * h;
    sys.evaluateResidual();
    for (size_t row = 0; row < n; ++row)
    {
      const double fd    = (plus[row] - sys.getResidual().getData()[row]) / (2 * h);
      const double error = std::abs(fd - entries[{row, column}]) / (1 + std::abs(fd));
      maximum            = std::max(maximum, error);
    }
    y[column]  = value;
    yp[column] = derivative;
  }
  sys.evaluateResidual();
  if (maximum > 1e-4)
    throw std::runtime_error("Assembled Jacobian finite-difference check failed: " + std::to_string(maximum));
  return maximum;
}

int main(int argc, char** argv)
{
  if (argc != 2 && argc != 3)
    throw std::invalid_argument("Usage: paraemt_9bus path/to/solver.json [trip]");
  const bool trip = argc == 3 && std::string(argv[2]) == "trip";
  if (argc == 3 && !trip)
    throw std::invalid_argument("Unknown event");
  const auto begin = std::chrono::steady_clock::now();
  auto       study = GridKit::EMT::parseStudyData(argv[1]);
  GridKit::EMT::configureCommonMath<double>(study);
  NineBus sys(study, std::filesystem::path(argv[1]).parent_path() / "network.state.json", trip);
  if (!sys.hasJacobian())
    throw std::runtime_error("Sparse Enzyme GridKit build required");
  sys.allocate();
  Ida<double, size_t> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setMaxSteps(study.max_steps);
  ida.setConsistentICType(study.consistent_ic_type);
  ida.configureSimulation();
  const double                               jacobian_initial = checkJacobian(sys, 0.0);
  GridKit::EMT::StateMonitor<double, size_t> state(sys, study);
  double                                     current_mismatch = 0;
  auto                                       record           = [&](double time)
  {
    state.write(time);
    current_mismatch = std::max(current_mismatch, sys.generatorCurrentMismatch());
  };
  std::cout << "DAE variables: " << sys.size() << ", Jacobian nonzeros: " << sys.getCsrJacobian()->getNnz() << std::endl;
  ida.initializeSimulation(0.0, false);
  sys.evaluateResidual();
  double initial_residual = 0;
  for (size_t i = 0; i < sys.size(); ++i)
    initial_residual = std::max(initial_residual, std::abs(sys.getResidual().getData()[i]));
  std::cout << "Initial maximum absolute residual (mixed physical units): " << initial_residual << std::endl;
  record(0);
  ida.runSimulation(1.0, study.dt_monitor, record);
  auto                      total       = ida.getStats();
  auto*                     pref        = sys.getSignal("pref_1");
  const double              pref_before = pref->read();
  const std::vector<double> before_event(sys.y().getData(), sys.y().getData() + sys.size());
  json                      projection = json::array();
  if (trip)
  {
    projection = sys.tripGenerator();
    std::cout << "Generator 1 ideal terminal opening at 1 s" << std::endl;
  }
  else
  {
    pref->init(pref_before - .02);
    std::cout << "Governor reference step: " << std::setprecision(15) << pref_before << " -> " << pref->read() << " at 1 s" << std::endl;
  }
  const std::vector<double> event_state(sys.y().getData(), sys.y().getData() + sys.size());
  if (trip)
    ida.getDefaultInitialCondition();
  // IDA recomputes algebraic values and derivatives, retaining the explicitly
  // supplied differential states (projected first for an ideal opening).
  ida.initializeSimulation(1.0);
  double event_state_change = 0;
  for (size_t i = 0; i < sys.size(); ++i)
    if (sys.tag()[i])
      event_state_change = std::max(event_state_change, std::abs(sys.y().getData()[i] - event_state[i]));
  if (event_state_change != 0)
    throw std::runtime_error("Event consistency solve changed a differential state");
  std::ofstream("event_state_limits.json") << json{{"time_s", 1.0},
                                                   {"before", before_event},
                                                   {"after_projection", event_state},
                                                   {"after_consistency", std::vector<double>(sys.y().getData(), sys.y().getData() + sys.size())},
                                                   {"differential", sys.tag()}}
                                                  .dump(2)
                                           << '\n';
  ida.runSimulation(study.tmax, study.dt_monitor, record);
  total                       += ida.getStats();
  const double jacobian_final  = checkJacobian(sys, study.tmax);
  sys.stopMonitor();
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
  std::ofstream("run.json") << json{{"simulator", "GridKit EMT"}, {"duration_s", study.tmax}, {"output_interval_s", study.dt_monitor}, {"rel_tol", study.rel_tol}, {"abs_tol", study.abs_tol}, {"mu", study.mu}, {"dae_variables", sys.size()}, {"jacobian_check_max_scaled_difference", {{"initial", jacobian_initial}, {"final", jacobian_final}}}, {"initial_max_absolute_residual_mixed_units", initial_residual}, {"max_generator_kcl_mismatch_A", current_mismatch}, {"wall_s", seconds}, {"event_differential_state_change", event_state_change}, {"event", {{"time_s", 1.0}, {"type", trip ? "ideal-generator-terminal-opening" : "governor-reference-step"}, {"generator_bus", 1}, {"pref_before", pref_before}, {"increment_pu", trip ? 0.0 : -.02}, {"explicit_state_projection", projection}}}, {"solver", {{"steps", total.num_steps_}, {"residual_evaluations", total.num_residual_evals_}, {"linear_setups", total.num_linear_decompositions_}, {"error_test_failures", total.num_error_test_fails_}, {"nonlinear_iterations", total.num_nonlinear_iters_}, {"nonlinear_convergence_failures", total.num_nonlinear_convergence_fails_}}}}.dump(2) << '\n';
  std::cout << "IDA steps: " << total.num_steps_ << ", residual evaluations: " << total.num_residual_evals_ << "\nComplete in " << seconds << " wall seconds\n";
}
