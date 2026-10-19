#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp>

#include "LineCircuit.hpp"

using namespace GridKit::Examples::EMT;
using System = GridKit::PhasorDynamics::SystemModel<double, size_t>;
using Solver = AnalysisManager::Sundials::Ida<double, size_t>;

// Save self-contained coefficients with each run for reproducible validation.
Json resolve(Json value, const std::filesystem::path& directory, bool line = false)
{
  auto base = directory;
  if (value.is_string())
  {
    const auto path = base / value.get<std::string>();
    value           = readJson(path);
    base            = path.parent_path();
  }
  if (line)
  {
    value["Yc"] = resolve(value.at("Yc"), base);
    auto hbase  = base;
    Json h      = value.at("H");
    if (h.is_string())
    {
      const auto path = base / h.get<std::string>();
      h               = readJson(path);
      hbase           = path.parent_path();
    }
    if (h.contains("fit"))
    {
      h["fit"] = resolve(h.at("fit"), hbase);
    }
    else
    {
      h["input"]  = resolve(h.at("input"), hbase);
      h["output"] = resolve(h.at("output"), hbase);
    }
    value["H"] = std::move(h);
  }
  return value;
}

void run(Json input, const std::filesystem::path& directory, const std::filesystem::path& destination)
{
  input["line"]          = resolve(input.at("line"), directory, true);
  const auto        data = CircuitData::parse(input, directory);
  const std::string name = input.value("name", "simulation");
  if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos)
  {
    throw std::invalid_argument("Run name must contain only letters, digits, '-' or '_'");
  }
  const auto output = destination / name;
  std::filesystem::create_directories(output);
  // A failed rerun must not retain a previous success marker or reference.
  std::filesystem::remove(output / "statistics.json");
  std::filesystem::remove(output / "reference.csv");
  std::ofstream(output / "input.json") << input.dump(2) << '\n';

  System     system(signalData(data.line));
  DelayPorts ports;
  for (size_t i = 0; i < 2 * data.line.K; ++i)
  {
    ports.source.push_back(system.getSignalNode(i));
    ports.delayed.push_back(system.getSignalNode(2 * data.line.K + i));
  }
  auto* circuit = new LineCircuit(data, std::move(ports));
  system.addComponent(circuit);
  const Complex jw{0, data.omega()};
  for (size_t i = 0; i < 2 * data.line.K; ++i)
  {
    if (data.line.delays[i % data.line.K] == 0)
    {
      continue;
    }
    const Complex amplitude = circuit->modalAmplitude(i);
    system.setSignalPrehistory(i, [amplitude, jw](double t)
                               { return std::real(amplitude * std::exp(jw * t)); },
                               [amplitude, jw](double t)
                               { return std::real(jw * amplitude * std::exp(jw * t)); });
  }
  if (!system.hasJacobian())
  {
    throw std::runtime_error("This example requires GridKit's Enzyme-enabled system Jacobian assembly (the example itself uses analytic derivatives)");
  }
  system.updateTime(data.t0, 0);
  system.allocate();
  system.stopMonitor();
  Solver solver(&system);
  if (system.hasSignalDelays() && !Solver::supportsStepControl())
  {
    throw std::runtime_error("LineDistributed requires IDAS step-check support");
  }
  // Reserve most of the requested error budget for history reconstruction.
  // In particular, delayed algebraic variables must resolve history_abs_tol.
  solver.setTolerance(data.rtol * .1);
  solver.setMaxSteps(1000000);
  solver.setMaxOrder(data.max_order);
  solver.setSuppressAlgebraicErrors(data.suppress_algebraic);
  GridKit::PhasorDynamics::configureSignalDelays(system, solver, {.rel_tol = data.rtol, .abs_tol = data.history_atol, .compress = data.compress});

  std::ofstream samples(output / "samples.csv"), steps(output / "steps.csv"), errors(output / "errors.csv");
  samples << std::setprecision(17);
  steps << std::setprecision(17) << "t,h,h_over_tau,history_error,segments,order\n";
  // Diagnostic: the variables dominating IDAS's weighted local error estimate.
  errors << std::setprecision(6) << "t,h,order,wrms,name1,w1,name2,w2,name3,w3,name4,w4\n";
  std::vector<std::string> names(system.size(), "?");
  for (size_t i = 0; i < circuit->localCount(); ++i)
  {
    names.at(circuit->globalIndex(i)) = circuit->variableName(i);
  }
  for (size_t i = 0; i < 2 * data.line.K; ++i)
  {
    if (data.line.delays[i % data.line.K] > 0)
    {
      names.at(system.getSignalNode(2 * data.line.K + i)->getVariableIndex()) = "delayed_" + std::to_string(i / data.line.K + 1) + "_k" + std::to_string(i % data.line.K + 1);
    }
  }
  const auto columns = circuit->columns();
  for (size_t i = 0; i < columns.size(); ++i)
  {
    samples << (i ? "," : "") << columns[i];
  }
  samples << '\n';
  auto observe = [&](double t)
  {
    const auto row = circuit->sample(t);
    for (size_t i = 0; i < row.size(); ++i)
    {
      samples << (i ? "," : "") << row[i];
    }
    samples << '\n';
  };
  double candidate_begin = data.t0, candidate_end = data.t0, max_step = 0;
  size_t max_segments = 0, accepted = 0;
  double min_delay = std::numeric_limits<double>::infinity();
  for (double delay : data.line.delays)
  {
    if (delay > 0)
    {
      min_delay = std::min(min_delay, delay);
    }
  }
  if (system.hasSignalDelays())
  {
    auto control   = solver.getStepControl();
    control.assess = [&, assess = control.assess](double begin, double end, const Solver::DenseOutput& dense)
    {
      candidate_begin = begin;
      candidate_end   = end;
      return assess(begin, end, dense);
    };
    control.finish = [&, finish = control.finish]()
    {
      finish();
      size_t segments = 0;
      for (size_t i = 0; i < 2 * data.line.K; ++i)
      {
        if (data.line.delays[i % data.line.K] > 0)
        {
          segments += system.signalHistory(i).segmentCount();
        }
      }
      const double h = candidate_end - candidate_begin;
      max_step       = std::max(max_step, h);
      max_segments   = std::max(max_segments, segments);
      ++accepted;
      steps << candidate_end << ',' << h << ',' << h / min_delay << ','
            << solver.getLastStepAssessment().error << ',' << segments << ',' << solver.getLastOrder() << '\n';
      if (accepted % 50 == 0)
      {
        const auto          ele = solver.getEstimatedLocalErrors();
        const auto          ewt = solver.getErrorWeights();
        std::vector<double> weighted(ele.size());
        double              sum = 0;
        for (size_t i = 0; i < ele.size(); ++i)
        {
          weighted[i] = std::abs(ele[i]) * ewt[i];
          sum += weighted[i] * weighted[i];
        }
        std::vector<size_t> order(ele.size());
        std::iota(order.begin(), order.end(), 0);
        std::partial_sort(order.begin(), order.begin() + 4, order.end(), [&](size_t a, size_t b)
                          { return weighted[a] > weighted[b]; });
        errors << candidate_end << ',' << h << ',' << solver.getLastOrder() << ',' << std::sqrt(sum / ele.size());
        for (size_t j = 0; j < 4; ++j)
        {
          errors << ',' << names[order[j]] << ',' << weighted[order[j]];
        }
        errors << '\n';
      }
    };
    solver.setStepControl(std::move(control));
  }
  solver.configureSimulation();
  solver.initializeSimulation(data.t0, false); // Exact zero or sinusoidal initial state and prehistory.
  observe(data.t0);

  // Land at event boundaries. Narrow transitions receive a local resolution cap;
  // smooth intervals retain the user's independent h_max (zero means uncapped).
  std::set<double> stops{data.t0, data.tf};
  for (const auto& event : data.events)
  {
    for (double t : {event.on, event.on + event.rise, event.off, event.off + event.fall})
    {
      if (t > data.t0 && t < data.tf)
      {
        stops.insert(t);
      }
    }
  }
  const auto start = std::chrono::steady_clock::now();
  double     begin = data.t0;
  for (double end : stops)
  {
    if (end <= begin)
    {
      continue;
    }
    double       cap    = data.hmax;
    const double middle = (begin + end) / 2;
    for (const auto& event : data.events)
    {
      double width = 0;
      if (middle > event.on && middle < event.on + event.rise)
      {
        width = event.rise;
      }
      if (middle > event.off && middle < event.off + event.fall)
      {
        width = event.fall;
      }
      if (width > 0)
      {
        cap = cap > 0 ? std::min(cap, width / 8) : width / 8;
      }
    }
    solver.setMaxStep(cap);
    solver.runSimulation(end, data.dt, observe);
    begin = end;
  }
  const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  const auto   stats   = solver.getStats();
  Json         report  = {{"seconds_including_csv", elapsed}, {"unknowns", system.size()}, {"steps", stats.num_steps_}, {"history_accepted_steps", accepted}, {"history_rejections", solver.getRejectedStepCount()}, {"error_test_failures_including_history", stats.num_error_test_fails_}, {"nonlinear_failures", stats.num_nonlinear_convergence_fails_}, {"residual_evaluations", stats.num_residual_evals_}, {"factorizations", stats.num_linear_decompositions_}, {"max_step", max_step}, {"max_history_segments", max_segments}};
  std::ofstream(output / "statistics.json") << report.dump(2) << '\n';
  std::cout << name << ": " << system.size() << " unknowns, " << stats.num_steps_
            << " steps, Complete in " << elapsed << " s\n";
}

int main(int argc, char** argv)
{
  try
  {
    if (argc < 2 || argc > 3)
    {
      throw std::invalid_argument("Usage: line_distributed input.json [results-directory]");
    }
    const auto path        = std::filesystem::absolute(argv[1]);
    const auto destination = argc == 3 ? std::filesystem::absolute(argv[2]) : path.parent_path() / "results";
    const auto input       = readJson(path);
    if (input.contains("runs"))
    {
      fields(input, {"defaults", "runs"});
      std::set<std::string> names;
      for (const auto& overrides : input.at("runs"))
      {
        Json merged = input.at("defaults");
        merged.merge_patch(overrides);
        if (!names.insert(merged.at("name").get<std::string>()).second)
        {
          throw std::invalid_argument("Duplicate run name");
        }
        run(merged, path.parent_path(), destination);
      }
    }
    else
    {
      run(input, path.parent_path(), destination);
    }
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "LineDistributed: " << error.what() << '\n';
    return 1;
  }
}
