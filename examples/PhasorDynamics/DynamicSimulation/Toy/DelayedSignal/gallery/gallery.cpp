#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
using namespace GridKit;
using namespace GridKit::PhasorDynamics;
using System = SystemModel<double, size_t>;
using Solver = AnalysisManager::Sundials::Ida<double, size_t>;
using Node   = SignalNode<double, size_t>;

class Feedback : public PhasorDynamics::Component<double, size_t>
{
  Node&                   source_;
  Node&                   delayed_;
  double                  gamma_;
  static constexpr double lambda = -0.4;
  static constexpr double k      = 1.3;

public:
  Feedback(Node& source, Node& delayed, double tau)
    : source_(source), delayed_(delayed), gamma_(k * std::exp(-lambda * tau) - lambda)
  {
    size_ = 1;
    source_.setAssigned();
  }

  int setGridKitComponentID(size_t id) override
  {
    gridkit_component_id_ = id;
    return 0;
  }

  int allocate() override
  {
    if (!allocated_)
    {
      allocateVectors(1);
    }
    tag_.assign(1, true);
    variable_indices_.assign(1, 0);
    residual_indices_.assign(1, 0);
    source_.link(y_.getData(), &getVariableIndex(0));
    allocated_ = true;
    return 0;
  }

  int verify() const override
  {
    return 0;
  }

  int initialize() override
  {
    y_.getData()[0]  = 1;
    yp_.getData()[0] = lambda;
    return 0;
  }

  int tagDifferentiable() override
  {
    tag_[0] = true;
    return 0;
  }

  int setAbsoluteTolerance(double tolerance) override
  {
    abs_tol_.getData()[0] = tolerance;
    return 0;
  }

  int evaluateResidual() override
  {
    f_.getData()[0] = yp_.getData()[0] + gamma_ * y_.getData()[0] - k * delayed_.read();
    return 0;
  }

  int evaluateJacobian() override
  {
    if (!J_rows_buffer_)
    {
      nnz_           = 2;
      J_rows_buffer_ = new size_t[2];
      J_cols_buffer_ = new size_t[2];
      J_vals_buffer_ = new double[2];
    }
    J_rows_buffer_[0] = J_rows_buffer_[1] = getResidualIndex(0);
    J_cols_buffer_[0]                     = getVariableIndex(0);
    J_cols_buffer_[1]                     = delayed_.getVariableIndex();
    J_vals_buffer_[0]                     = gamma_ + alpha_;
    J_vals_buffer_[1]                     = -k;
    return constructCoo();
  }
};

// A manufactured algebraic source; delay and integration code are linked GridKit.
class Source : public Component<double, size_t>
{
  Node&                         node_;
  std::function<double(double)> function_;

public:
  Source(Node& n, std::function<double(double)> f)
    : node_(n), function_(f)
  {
    size_ = 1;
    node_.setAssigned();
  }

  int setGridKitComponentID(size_t id) override
  {
    gridkit_component_id_ = id;
    return 0;
  }

  int allocate() override
  {
    if (!allocated_)
      allocateVectors(1);
    tag_.assign(1, false);
    variable_indices_.assign(1, 0);
    residual_indices_.assign(1, 0);
    node_.link(y_.getData(), &getVariableIndex(0));
    allocated_ = true;
    return 0;
  }

  int verify() const override
  {
    return 0;
  }

  int initialize() override
  {
    y_.getData()[0]  = function_(0);
    yp_.getData()[0] = 0;
    return 0;
  }

  int tagDifferentiable() override
  {
    tag_[0] = false;
    return 0;
  }

  int setAbsoluteTolerance(double x) override
  {
    abs_tol_.getData()[0] = x;
    return 0;
  }

  int evaluateResidual() override
  {
    f_.getData()[0] = y_.getData()[0] - function_(time_);
    return 0;
  }

  int evaluateJacobian() override
  {
    if (!J_rows_buffer_)
    {
      nnz_           = 1;
      J_rows_buffer_ = new size_t[1];
      J_cols_buffer_ = new size_t[1];
      J_vals_buffer_ = new double[1];
    }
    J_rows_buffer_[0] = getResidualIndex(0);
    J_cols_buffer_[0] = getVariableIndex(0);
    J_vals_buffer_[0] = 1;
    return constructCoo();
  }
};

int main(int argc, char** argv)
{
  try
  {
    if (argc != 11)
      throw std::runtime_error("name kind tau tf rtol compress fixed_h width dt history_rtol");
    const std::string             name = argv[1], kind = argv[2];
    const double                  tau = std::stod(argv[3]), tf = std::stod(argv[4]), tol = std::stod(argv[5]);
    const bool                    compress = std::stoi(argv[6]);
    const double                  h = std::stod(argv[7]), width = std::stod(argv[8]), dt = std::stod(argv[9]), htol = std::stod(argv[10]);
    std::function<double(double)> exact = [&](double t)
    {
      if (kind == "feedback")
        return std::exp(-0.4 * t);
      if (kind == "pulse")
        return 1 + 0.6 * std::exp(-std::pow((t - 2) / width, 2));
      if (kind == "chirp")
        return std::cos(t + 0.4 * t * t);
      return std::cos(t);
    };
    SystemModelData<> data;
    data.signal = {{.name = "u", .signal_id = 0}, {.name = "delayed", .signal_id = 1, .source = 0, .delay = tau}};
    System system(data);
    if (kind == "feedback")
      system.addComponent(new Feedback(*system.getSignalNode(0), *system.getSignalNode(1), tau));
    else
      system.addComponent(new Source(*system.getSignalNode(0), exact));
    system.signalHistory(0).setPrehistory(exact);
    system.allocate();
    system.stopMonitor();
    system.setSignalHistoryTolerance(htol, htol * .01, compress);
    Solver solver(&system);
    solver.setTolerance(tol, tol * .01);
    solver.setMaxSteps(1000000);
    if (h > 0)
      solver.setFixedStep(h);
    std::ofstream steps(name + "_steps.csv"), samples(name + "_samples.csv"), attempts(name + "_attempts.csv"), midpoints(name + "_midpoints.csv");
    for (auto* s : {&steps, &samples, &attempts, &midpoints})
      *s << std::setprecision(17);
    steps << "t,h,h_over_tau,segments,u,z,error_estimate\n";
    samples << "t,u,z,exact_u,exact_z\n";
    attempts << "t,h,error_estimate,rejected\n";
    midpoints << "t,u,z,exact_u,exact_z\n";
    double begin = 0, end = 0, error = 0, mid_u = 0, mid_z = 0;
    solver.setStepControl({[&](double t, bool consistent)
                           { system.initializeSignalHistory(t, consistent); },
                           [&](double a, double b, const Solver::DenseOutput& dense)
                           {begin=a;end=b;error=system.assessSignalHistory(a,b,dense);attempts<<end<<','<<end-begin<<','<<error<<','<<(error>1)<<'\n';dense((begin+end)/2);mid_u=system.getSignalNode(0)->read();mid_z=system.getSignalNode(1)->read();return error; },
                           [&]()
                           { system.acceptSignalHistory(); },
                           [&]()
                           {system.finishSignalHistory();const double tm=(begin+end)/2;midpoints<<tm<<','<<mid_u<<','<<mid_z<<','<<exact(tm)<<','<<exact(tm-tau)<<'\n';steps<<end<<','<<end-begin<<','<<(end-begin)/tau<<','<<system.signalHistory(0).segmentCount()<<','<<system.getSignalNode(0)->read()<<','<<system.getSignalNode(1)->read()<<','<<error<<'\n'; },
                           [&]()
                           { system.discardSignalHistory(); }});
    solver.configureSimulation();
    solver.initializeSimulation(0);
    auto observe = [&](double t)
    { samples << t << ',' << system.getSignalNode(0)->read() << ',' << system.getSignalNode(1)->read() << ',' << exact(t) << ',' << exact(t - tau) << '\n'; };
    observe(0);
    solver.runSimulation(tf, dt, observe);
    std::cout << name << " rejected=" << solver.getRejectedStepCount() << '\n';
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
