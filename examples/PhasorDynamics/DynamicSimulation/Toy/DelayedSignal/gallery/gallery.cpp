// Manufactured smooth sources, integrated by the actual GridKit delay solver.
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>

using namespace GridKit;
using namespace GridKit::PhasorDynamics;
using System = SystemModel<double, size_t>;
using Solver = AnalysisManager::Sundials::Ida<double, size_t>;
using Node   = SignalNode<double, size_t>;

// Like FunctionSignalSource, this prescribed input is evaluated at the queried time.
class Source : public Component<double, size_t>
{
  Node&                         node_;
  std::function<double(double)> function_;
  double                        value_{0};
  size_t                        index_{INVALID_INDEX<size_t>};

public:
  Source(Node& node, std::function<double(double)> function)
    : node_(node), function_(std::move(function))
  {
    size_ = 0;
    node_.setAssigned();
  }

  int setGridKitComponentID(size_t id) override
  {
    gridkit_component_id_ = id;
    return 0;
  }

  int allocate() override
  {
    node_.link(&value_, &index_);
    allocated_ = true;
    return 0;
  }

  int verify() const override
  {
    return 0;
  }

  int initialize() override
  {
    value_ = function_(time_);
    return 0;
  }

  int tagDifferentiable() override
  {
    return 0;
  }

  int setAbsoluteTolerance(double) override
  {
    return 0;
  }

  int evaluateResidual() override
  {
    return 0;
  }

  int evaluateJacobian() override
  {
    return constructCoo();
  }

  void updateTime(double t, double alpha) override
  {
    Component<double, size_t>::updateTime(t, alpha);
    value_ = function_(t);
  }
};

struct Waveform
{
  std::function<double(double)> value;
  std::function<double(double)> derivative;
  std::vector<double>           corners;
};

Waveform inputFunction(const nlohmann::json& run)
{
  const auto   shape  = run.at("shape").get<std::string>();
  const double width  = run.at("width");
  const double center = run.value("center", 2.0);
  auto         pulse  = [](double t, double c, double w)
  { return std::exp(-std::pow((t - c) / w, 2)); };
  auto slope = [pulse](double t, double c, double w)
  { return -2 * (t - c) * pulse(t, c, w) / (w * w); };
  if (shape == "pulse")
  {
    return {[=](double t)
            { return 1 + 0.6 * pulse(t, center, width); },
            [=](double t)
            { return 0.6 * slope(t, center, width); },
            {}};
  }
  if (shape == "ramp" || shape == "rise" || shape == "fall")
  {
    const double scale = shape == "ramp" ? 1 / width : std::log(9.0) / width;
    const double sign  = shape == "fall" ? -1 : 1;
    return {[=](double t)
            { return 1 + 0.3 * (1 + sign * std::tanh(scale * (t - center))); },
            [=](double t)
            { const double y = std::tanh(scale * (t - center)); return 0.3 * sign * scale * (1 - y * y); },
            {}};
  }
  if (shape == "asymmetric")
  {
    const double fall_width = run.at("fall_width");
    const double fall_time  = run.at("fall_time");
    return {[=](double t)
            { return 1 + 0.3 * (std::tanh(std::log(9.0) * (t - center) / width) - std::tanh(std::log(9.0) * (t - fall_time) / fall_width)); },
            [=](double t)
            { const double rise = std::tanh(std::log(9.0) * (t - center) / width);
                             const double fall = std::tanh(std::log(9.0) * (t - fall_time) / fall_width);
                             return 0.3 * std::log(9.0) * ((1 - rise * rise) / width - (1 - fall * fall) / fall_width); },
            {}};
  }
  if (shape == "features")
  {
    return {[=](double t)
            { return 1 + 0.25 * pulse(t, 2, 0.3) + 0.6 * pulse(t, 4, width); },
            [=](double t)
            { return 0.25 * slope(t, 2, 0.3) + 0.6 * slope(t, 4, width); },
            {}};
  }
  if (shape == "packet")
  {
    return {[=](double t)
            { return 1 + 0.5 * pulse(t, 3, width) * std::sin(20 * (t - 3)); },
            [=](double t)
            { return 0.5 * slope(t, 3, width) * std::sin(20 * (t - 3))
                     + 10 * pulse(t, 3, width) * std::cos(20 * (t - 3)); },
            {}};
  }
  if (shape == "sine" || shape == "chirp")
  {
    const double frequency = run.at("frequency");
    const double rate      = run.value("chirp_rate", 0.0);
    const double pi        = std::numbers::pi;
    return {[=](double t)
            { return std::sin(2 * pi * (frequency * t + 0.5 * rate * t * t)); },
            [=](double t)
            { return 2 * pi * (frequency + rate * t)
                     * std::cos(2 * pi * (frequency * t + 0.5 * rate * t * t)); },
            {}};
  }
  if (shape == "triangle")
  {
    return {[](double t)
            { if (t < 0.5) return 0.0;
                            if (t < 1.0) return 2 * (t - 0.5);
                            if (t < 1.5) return 2 * (1.5 - t);
                            return 0.0; },
            [](double t)
            { if (t < 0.5 || t >= 1.5) return 0.0;
                            return t < 1.0 ? 2.0 : -2.0; },
            {0.5, 1.0, 1.5}};
  }
  if (shape == "corner")
  {
    return {[](double t)
            { return t < 0.5 ? 0.0 : (t < 1.0 ? t - 0.5 : 0.5 + 0.2 * (t - 1.0)); },
            [](double t)
            { return t < 0.5 ? 0.0 : (t < 1.0 ? 1.0 : 0.2); },
            {0.5, 1.0}};
  }
  throw std::invalid_argument("Unknown input shape: " + shape);
}

int main(int argc, char** argv)
{
  try
  {
    if (argc != 3)
      throw std::invalid_argument("Usage: gallery run.json output_prefix");
    std::ifstream     input(argv[1]);
    const auto        run  = nlohmann::json::parse(input);
    const std::string name = argv[2];
    const double      tau = run.at("tau"), tolerance = run.at("rel_tol");
    const auto        waveform = inputFunction(run);
    const auto&       exact    = waveform.value;

    SystemModelData<> data;
    data.signal = {{.name = "input", .signal_id = 0},
                   {.name = "delayed", .signal_id = 1, .source = 0, .delay = tau}};
    System system(data);
    system.addComponent(new Source(*system.getSignalNode(0), exact));
    system.setSignalPrehistory(0, exact, waveform.derivative);
    system.allocate();
    system.stopMonitor();
    Solver solver(&system);
    solver.setTolerance(tolerance, tolerance * 0.01);
    solver.setMaxSteps(1000000);
    solver.setMinStep(run.at("h_min"));
    solver.setMaxStep(run.value("h_max", 0.0));
    configureSignalDelays(system, solver, {.rel_tol = tolerance, .abs_tol = tolerance * 0.01, .compress = run.at("compress").get<bool>()});

    // Compare storage on identical accepted knots, without affecting the solve.
    // Index 0 is uncompressed; index 1 enables the actual history compressor.
    std::array<SignalHistory<double>, 2> histories;
    for (size_t i = 0; i < histories.size(); ++i)
    {
      histories[i].setOptions({.rel_tol = tolerance, .abs_tol = tolerance * 0.01, .compress = i == 1});
      histories[i].setPrehistory(exact, waveform.derivative);
    }
    const std::array<double, 1> delays{tau};

    std::ofstream steps(name + "_steps.csv"), samples(name + "_samples.csv"),
        attempts(name + "_attempts.csv"), midpoints(name + "_midpoints.csv"), segments(name + "_segments.csv");
    for (auto* stream : {&steps, &samples, &attempts, &midpoints, &segments})
      *stream << std::setprecision(17);
    steps << "t,h,h_over_tau,segments,u,z,error_estimate,exact_z,native_error_failures,nonlinear_failures,uncompressed_segments,compressed_segments\n";
    samples << "t,u,z,exact_u,exact_z\n";
    midpoints << "t,u,z,exact_u,exact_z\n";
    attempts << "t,h,error_estimate,rejected,check,error_order\n";
    segments << "begin,end,c0,c1,c2,c3\n";
    long   native_before = 0, nonlinear_before = 0;
    size_t rejected_total = 0, rejected_before = 0;

    // Decorate the installed control only to record diagnostics.
    auto   control = solver.getStepControl();
    double begin = 0, end = 0, error = 0, mid_u = 0, mid_z = 0;
    control.assess = [&, assess = control.assess](double a, double b, const Solver::DenseOutput& dense)
    {
      begin                  = a;
      end                    = b;
      const auto assessment  = assess(a, b, dense);
      error                  = assessment.error;
      const bool rejected    = !std::isfinite(error) || error > 1;
      rejected_total        += rejected;
      const int check        = solver.getLastStepDescription().starts_with("source") ? 0 : 1;
      attempts << end << ',' << end - begin << ',' << error << ',' << rejected
               << ',' << check << ',' << assessment.error_order << '\n';
      dense((begin + end) / 2);
      mid_u = system.getSignalNode(0)->read();
      mid_z = system.getSignalNode(1)->read();
      return assessment;
    };
    control.accept = [&, accept = control.accept]()
    {
      accept();
      const auto& segment = system.signalHistory(0).segments().back();
      for (auto& history : histories)
      {
        history.prepare(segment.end, segment.coefficients[3]);
        history.accept();
      }
      segments << segment.begin << ',' << segment.end;
      for (double coefficient : segment.coefficients)
      {
        segments << ',' << coefficient;
      }
      segments << '\n';
    };
    control.finish = [&, finish = control.finish]()
    {
      finish();
      for (auto& history : histories)
      {
        history.maintain(delays);
      }
      if (histories[run.at("compress").get<bool>() ? 1 : 0].segmentCount()
          != system.signalHistory(0).segmentCount())
      {
        throw std::logic_error("Recorded history comparison differs from the solver history");
      }
      const double t = (begin + end) / 2;
      midpoints << t << ',' << mid_u << ',' << mid_z << ',' << exact(t) << ',' << exact(t - tau) << '\n';
      steps << end << ',' << end - begin << ',' << (end - begin) / tau << ','
            << system.signalHistory(0).segmentCount() << ','
            << system.getSignalNode(0)->read() << ',' << system.getSignalNode(1)->read()
            << ',' << error << ',' << exact(end - tau);
      const auto stats = solver.getStats();
      steps << ',' << native_before + stats.num_error_test_fails_ - static_cast<long>(rejected_total - rejected_before)
            << ',' << nonlinear_before + stats.num_nonlinear_convergence_fails_
            << ',' << histories[0].segmentCount() << ',' << histories[1].segmentCount() << '\n';
    };
    solver.setStepControl(std::move(control));
    solver.configureSimulation();
    solver.initializeSimulation(0);
    for (auto& history : histories)
    {
      history.initialize(0, system.getSignalNode(0)->read());
    }
    auto observe = [&](double t)
    {
      samples << t << ',' << system.getSignalNode(0)->read() << ','
              << system.getSignalNode(1)->read() << ',' << exact(t) << ',' << exact(t - tau) << '\n';
    };
    observe(0);
    auto complete_phase = [&]()
    {
      const auto stats  = solver.getStats();
      native_before    += stats.num_error_test_fails_ - static_cast<long>(rejected_total - rejected_before);
      nonlinear_before += stats.num_nonlinear_convergence_fails_;
      rejected_before   = rejected_total;
    };
    if (run.value("restart_corners", false))
    {
      for (double corner : waveform.corners)
      {
        solver.runSimulation(corner, run.at("dt_monitor"), observe);
        complete_phase();
        system.restartSignalHistory(0);
        for (auto& history : histories)
        {
          history.restartInterpolation();
        }
        solver.initializeSimulation(corner, false);
      }
    }
    solver.runSimulation(run.at("tf"), run.at("dt_monitor"), observe);
    complete_phase();
    std::cout << name << " rejected=" << rejected_total
              << " native_error_failures=" << native_before
              << " nonlinear_failures=" << nonlinear_before << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
