#include <cmath>
#include <limits>
#include <map>

#include <GridKit/Model/EMT/Component/Controller/DCLink/DcLink.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using Link   = EMT::Controller::DcLink<double, size_t>;
  using Signal = Link::SignalT;

  template <typename Action>
  bool rejects(Action action)
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

  Testing::TestOutcome validation()
  {
    Testing::TestStatus success = true;
    Link::ModelDataT    data;
    success *= rejects([&]
                       { Link link(data); });
    for (const double value : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      data.parameters[EMT::Controller::DcLinkParameters::C]  = value;
      success                                               *= rejects([&]
                         { Link link(data); });
    }
    data.parameters[EMT::Controller::DcLinkParameters::C]  = true;
    success                                               *= rejects([&]
                       { Link link(data); });
    data.parameters[EMT::Controller::DcLinkParameters::C]  = 0.02;
    Link link(data);
    success *= rejects([&]
                       { link.allocate(); });
    Signal source, sink, alias;
    source.bindConstant(80.0);
    sink.bindConstant(40.0);
    link.attachInput(&source, &sink);
    link.assignOutput(Link::Outputs::vdc, &alias);
    link.allocate();
    success *= link.verify() == 0;
    link.initialize({{Link::Outputs::vdc, 600.0}});
    success *= alias.read() == 600.0 && source.read() == 80.0 && sink.read() == 40.0;
    success *= rejects([&]
                       { link.attachInput(&source, &sink); });
    success *= rejects([&]
                       { link.initializeState({{"vdc", std::numeric_limits<double>::infinity()}}); });
    success *= rejects([&]
                       { link.initializeState({{"typo", 1.0}}); });
    success *= link.initializationPorts().inputs.empty();
    return success.report("DCLink parameter, wiring and initialization contracts");
  }

  Testing::TestOutcome derivatives()
  {
    Testing::TestStatus success = true;
    Link::ModelDataT    data;
    data.parameters[EMT::Controller::DcLinkParameters::C] = 0.02;
    Link   link(data);
    Signal source, sink;
    double external = 3.0;
    size_t column   = 7;
    source.set(&external, &column);
    sink.setComputed([&]
                     { return 0.1 * link.outputSignal(Link::Outputs::vdc).read() + external * external; },
                     [&](Signal::GradientT& gradient, double scale)
                     {
                       link.outputSignal(Link::Outputs::vdc).appendGradient(gradient, 0.1 * scale);
                       source.appendGradient(gradient, 2.0 * external * scale);
                     });
    link.attachInput(&source, &sink);
    link.allocate();
    link.assignGlobalIndices(2);
    link.initialize({{Link::Outputs::vdc, 600.0}});
    link.yp().getData()[0] = -3300.0;
    auto residual          = [&]
    { link.evaluateResidual(); return link.getResidual().getData()[0]; };
    success *= std::abs(residual()) < 1e-12;
    for (const auto& scales : {std::pair{1.0, 0.0}, std::pair{0.0, 1.0}, std::pair{2.0, 3.0}, std::pair{0.0, 0.0}, std::pair{1.0, 0.0}})
    {
      std::map<size_t, double> jacobian;
      for (const auto& entry : link.jacobianEntries(scales.first, scales.second))
        jacobian[entry.column] += entry.value;
      constexpr double h          = 1e-3;
      auto             difference = [&](double& value, double scale)
      {
        value            += h * scale;
        const auto plus   = residual();
        value            -= 2 * h * scale;
        const auto minus  = residual();
        value            += h * scale;
        return (plus - minus) / (2 * h);
      };
      const auto owned  = difference(link.y().getData()[0], scales.first) + difference(link.yp().getData()[0], scales.second);
      success          *= std::abs(jacobian[2] - owned) < 1e-9;
      success          *= std::abs(jacobian[7] - difference(external, scales.first)) < 1e-9;
    }
    const auto mass  = link.jacobianEntries(0.0, 1.0);
    success         *= mass.size() == 1 && mass[0].column == 2 && mass[0].value == -0.02;
    return success.report("DCLink residual and composed Jacobian against finite differences");
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  results += validation();
  results += derivatives();
  return results.summary();
}
