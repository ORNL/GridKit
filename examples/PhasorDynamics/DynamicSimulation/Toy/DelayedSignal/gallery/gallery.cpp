// Manufactured smooth sources, integrated by the actual GridKit delay solver.
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

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

std::function<double(double)> inputFunction(const nlohmann::json& run)
{
  const auto   shape = run.at("shape").get<std::string>();
  const double width = run.at("width");
  auto         pulse = [](double t, double center, double duration)
  { return std::exp(-std::pow((t - center) / duration, 2)); };
  if (shape == "pulse")
    return [=](double t)
    { return 1 + 0.6 * pulse(t, 2, width); };
  if (shape == "ramp")
    return [=](double t)
    { return 1 + 0.3 * (1 + std::tanh((t - 2) / width)); };
  if (shape == "features")
    return [=](double t)
    { return 1 + 0.25 * pulse(t, 2, 0.3) + 0.6 * pulse(t, 4, width); };
  if (shape == "packet")
    return [=](double t)
    { return 1 + 0.5 * pulse(t, 3, width) * std::sin(20 * (t - 3)); };
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
    const auto        exact = inputFunction(run);

    SystemModelData<> data;
    data.signal = {{.name = "input", .signal_id = 0},
                   {.name = "delayed", .signal_id = 1, .source = 0, .delay = tau}};
    System system(data);
    system.addComponent(new Source(*system.getSignalNode(0), exact));
    system.setSignalPrehistory(0, exact);
    system.allocate();
    system.stopMonitor();
    Solver solver(&system);
    solver.setTolerance(tolerance, tolerance * 0.01);
    solver.setMaxSteps(1000000);
    configureSignalDelays(system, solver, {.rel_tol = tolerance, .abs_tol = tolerance * 0.01, .compress = run.at("compress").get<bool>()});

    std::ofstream steps(name + "_steps.csv"), samples(name + "_samples.csv"),
        attempts(name + "_attempts.csv"), midpoints(name + "_midpoints.csv");
    for (auto* stream : {&steps, &samples, &attempts, &midpoints})
      *stream << std::setprecision(17);
    steps << "t,h,h_over_tau,segments,u,z,error_estimate,exact_z\n";
    samples << "t,u,z,exact_u,exact_z\n";
    midpoints << "t,u,z,exact_u,exact_z\n";
    attempts << "t,h,error_estimate,rejected\n";

    // Decorate the installed control only to record diagnostics.
    auto   control = solver.getStepControl();
    double begin = 0, end = 0, error = 0, mid_u = 0, mid_z = 0;
    control.assess = [&, assess = control.assess](double a, double b, const Solver::DenseOutput& dense)
    {
      begin = a;
      end   = b;
      error = assess(a, b, dense);
      attempts << end << ',' << end - begin << ',' << error << ',' << (error > 1) << '\n';
      dense((begin + end) / 2);
      mid_u = system.getSignalNode(0)->read();
      mid_z = system.getSignalNode(1)->read();
      return error;
    };
    control.finish = [&, finish = control.finish]()
    {
      finish();
      const double t = (begin + end) / 2;
      midpoints << t << ',' << mid_u << ',' << mid_z << ',' << exact(t) << ',' << exact(t - tau) << '\n';
      steps << end << ',' << end - begin << ',' << (end - begin) / tau << ','
            << system.signalHistory(0).segmentCount() << ','
            << system.getSignalNode(0)->read() << ',' << system.getSignalNode(1)->read()
            << ',' << error << ',' << exact(end - tau) << '\n';
    };
    solver.setStepControl(std::move(control));
    solver.configureSimulation();
    solver.initializeSimulation(0);
    auto observe = [&](double t)
    {
      samples << t << ',' << system.getSignalNode(0)->read() << ','
              << system.getSignalNode(1)->read() << ',' << exact(t) << ',' << exact(t - tau) << '\n';
    };
    observe(0);
    solver.runSimulation(run.at("tf"), run.at("dt_monitor"), observe);
    std::cout << name << " rejected=" << solver.getRejectedStepCount() << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
