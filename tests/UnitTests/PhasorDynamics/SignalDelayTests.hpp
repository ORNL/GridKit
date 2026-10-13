#pragma once

#include <cmath>
#include <limits>

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelay.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  class SignalDelayTests
  {
    using History = PhasorDynamics::SignalHistory<double>;

    class InitializationWriter : public PhasorDynamics::Component<double, size_t>
    {
      PhasorDynamics::SignalNode<double, size_t>& input_;

    public:
      explicit InitializationWriter(PhasorDynamics::SignalNode<double, size_t>& input)
        : input_(input)
      {
      }

      int setGridKitComponentID(size_t id) override
      {
        gridkit_component_id_ = id;
        return 0;
      }

      int allocate() override
      {
        return 0;
      }

      int verify() const override
      {
        return 0;
      }

      int initialize() override
      {
        input_.init(2.4);
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
    };

    template <class Function>
    static bool throws(Function f)
    {
      try
      {
        f();
      }
      catch (const std::exception&)
      {
        return true;
      }
      return false;
    }

  public:
    TestOutcome prehistoryDerivative()
    {
      TestStatus success = true;
      History    history;
      auto       value = [](double t)
      { return std::sin(20 * t); };
      auto derivative = [](double t)
      { return 20 * std::cos(20 * t); };
      history.setPrehistory(value);
      success *= std::abs(history.initialDerivative(-8) - derivative(-8)) < 1e-5;
      history.setPrehistory(value, derivative);
      success *= history.initialDerivative(-8) == derivative(-8);
      history.initialize(0, 0);
      success *= throws([&]
                        { history.setPrehistory(value, derivative); });
      history.reset();
      history.setPrehistory(value, [](double)
                            { return std::numeric_limits<double>::quiet_NaN(); });
      success *= throws([&]
                        { history.initialDerivative(-1); });
      success *= throws([&]
                        { history.setPrehistory({}, derivative); });
      return success.report(__func__);
    }

    TestOutcome derivativeCorner()
    {
      TestStatus success = true;
      History    history;
      history.setOptions({.rel_tol = 0, .abs_tol = 1e-8, .compress = true});
      history.setPrehistory([](double t)
                            { return t; },
                            [](double)
                            { return 1.0; });
      history.initialize(0, 0);
      history.prepare(1, 1);
      history.accept();
      history.restartInterpolation();
      // A continuous triangular peak changes slope from +1 to -1.
      const auto sample  = history.sample(1.2, 0.1);
      success           *= std::abs(sample.constant + sample.weight * 0.8 - 0.9) < 1e-14;
      history.prepare(1.2, 0.8);
      history.discard();
      success *= history.time() == 1 && history.segmentCount() == 1;
      history.prepare(1.2, 0.8);
      history.accept();
      success *= history.segments().back().break_before;
      for (int i = 1; i <= 20; ++i)
      {
        const double t = 1.2 + i * 0.1;
        history.prepare(t, 2 - t);
        history.accept();
        history.maintain(std::array{4.0});
      }
      success              *= std::abs(history.sample(3.2, 2.3).constant - 0.9) < 1e-12;
      success              *= std::abs(history.sample(3.2, 2.1).constant - 0.9) < 1e-12;
      bool retained_corner  = false;
      for (const auto& segment : history.segments())
      {
        retained_corner |= segment.begin == 1 && segment.break_before;
        success         *= !(segment.begin < 1 && segment.end > 1);
      }
      success *= retained_corner;
      return success.report(__func__);
    }

    TestOutcome interpolation()
    {
      TestStatus success = true;
      for (int degree = 0; degree <= 3; ++degree)
      {
        auto polynomial = [degree](double t)
        { return std::pow(t, degree); };
        auto derivative = [degree](double t)
        { return degree == 0 ? 0 : degree * std::pow(t, degree - 1); };
        const std::array<History::Sample, 3> anchors{{{-0.7, polynomial(-0.7), derivative(-0.7)}, {-0.3, polynomial(-0.3), derivative(-0.3)}, {0.2, polynomial(0.2), derivative(0.2)}}};
        for (double q : {0.2, 0.31, 0.59, 0.8})
        {
          const auto a             = History::trialSample(anchors, 0.8, q);
          success                 *= std::abs(a.constant + a.weight * polynomial(0.8) - polynomial(q)) < 1e-12;
          const double epsilon     = 1e-6;
          const double derivative  = ((a.constant + a.weight * (polynomial(0.8) + epsilon))
                                     - (a.constant + a.weight * polynomial(0.8)))
                                    / epsilon;
          success *= std::abs(derivative - a.weight) < 1e-9;
        }
      }
      return success.report(__func__);
    }

    TestOutcome historyLifecycle()
    {
      TestStatus success = true;
      History    h;
      h.setPrehistory([](double t)
                      { return 2 + 3 * t; });
      h.initialize(0, 2);
      success *= h.sample(0, 0.2).constant == 1.4;
      h.prepare(1, 5);
      const auto   trial   = h.sample(1, 0.4);
      const double before  = trial.constant + trial.weight * 5;
      success             *= std::abs(before - 3.8) < 1e-10;
      h.discard();
      success *= h.time() == 0 && h.segmentCount() == 0;
      h.prepare(1, 5);
      h.accept();
      success *= std::abs(h.sample(1, 0.4).constant - before) < 1e-12;
      success *= h.sample(1, 0.4).weight == 0;
      success *= throws([&]
                        { h.begin(0.5); });
      h.begin(1);
      success *= h.time() == 1;
      h.prepare(2, 8);
      for (double delay : {0.25, 0.75})
      {
        std::array<double, 6> values;
        for (size_t i = 0; i < values.size(); ++i)
        {
          values[i] = 2 + 3 * (1 + History::dense_nodes[i] - delay);
        }
        success   *= h.estimateDelayed(1, 2, delay, values) < 1e-4;
        values[2] += 0.1;
        success   *= h.estimateDelayed(1, 2, delay, values) > 1;
      }
      h.reset();
      success *= !h.initialized() && h.hasPrehistory();
      return success.report(__func__);
    }

    TestOutcome compression()
    {
      TestStatus success = true;
      History    full, compact;
      full.setOptions({.rel_tol = 0, .abs_tol = 1e-4});
      compact.setOptions({.rel_tol = 0, .abs_tol = 1e-4, .compress = true});
      full.initialize(0, 0);
      compact.initialize(0, 0);
      success *= throws([&]
                        { compact.setOptions({}); });
      const std::array<double, 2> delays{0.001, 1.0};
      bool                        reduced = false;
      for (int i = 1; i <= 400; ++i)
      {
        const double t = i * 0.005;
        full.prepare(t, std::sin(t));
        compact.prepare(t, std::sin(t));
        if (i == 1)
        {
          const double dense  = 0.375 * std::sin(t) + 9e-5;
          success            *= full.estimate(t / 2, dense) < 1;
          success            *= compact.estimate(t / 2, dense) > 1;
        }
        full.accept();
        compact.accept();
        const double accepted = compact.sample(t, delays[0]).constant;
        full.maintain(delays);
        compact.maintain(delays);
        success *= compact.sample(t, delays[0]).constant == accepted;
        reduced |= compact.segmentCount() < full.segmentCount();
        for (int j = 1; j <= 40; ++j)
        {
          const double delay  = std::min(t, 1.0) * j / 40;
          success            *= std::abs(full.sample(t, delay).constant - compact.sample(t, delay).constant) <= 2.5e-5 + 1e-12;
        }
        for (const auto& segment : compact.segments())
        {
          success *= segment.error <= 2.5e-5;
        }
      }
      success *= reduced;
      success *= full.segmentCount() <= 202;
      success *= throws([&]
                        { compact.sample(2, 1.5); });
      return success.report(__func__);
    }

    TestOutcome continuityAndHiddenFeature()
    {
      TestStatus success = true;
      History    history;
      history.initialize(0, 1);
      for (int i = 1; i <= 200; ++i)
      {
        const double t = i / 200.0;
        history.prepare(t, 1 + std::exp(-std::pow((t - 0.47) / 0.01, 2)));
        history.accept();
      }
      const auto& segments = history.segments();
      for (size_t i = 1; i < segments.size(); ++i)
      {
        const auto&  a      = segments[i - 1];
        const auto&  b      = segments[i];
        const double left   = 3 * (a.coefficients[3] - a.coefficients[2]) / (a.end - a.begin);
        const double right  = 3 * (b.coefficients[1] - b.coefficients[0]) / (b.end - b.begin);
        success            *= std::abs(left - right) < 1e-10;
      }
      history.prepare(1.1, 1);
      // The pulse falls between the six dense-output nodes. Checking every
      // history polynomial must still reject a constant delayed trajectory.
      success *= history.estimateDelayed(1, 2, 1, {1, 1, 1, 1, 1, 1}) > 1;

      History constant;
      constant.setOptions({.rel_tol = 0, .abs_tol = 1e-7});
      constant.initialize(0, 1);
      constant.prepare(1, 1);
      constant.accept();
      constant.prepare(2, 1);
      std::array<double, 6> quintic;
      for (size_t i = 0; i < quintic.size(); ++i)
      {
        quintic[i] = 1 + 1e-4 * std::pow(2 * History::dense_nodes[i] - 1, 5);
      }
      // (2x-1)^5 has alternating +/-1 Bernstein coefficients. Check the
      // fifth-degree bound, including subdivision across accepted/trial history.
      for (double delay : {0.5, 1.0})
      {
        success *= std::abs(constant.estimateDelayed(1, 2, delay, quintic) - 1000) < 1e-5;
      }
      return success.report(__func__);
    }

    TestOutcome component()
    {
      using namespace PhasorDynamics;
      TestStatus                      success = true;
      SignalDelayData<double, size_t> data;
      data.parameters[SignalDelayParameters::tau] = 0.02;
      History                    history;
      SignalNode<double, size_t> source, output;
      double                     input        = 2.3;
      size_t                     source_index = 1;
      source.link(&input, &source_index);
      SignalDelay<double, size_t> delay(data, history);
      delay.getPorts().in.template port<SignalDelayInputs::u>().connect(&source);
      delay.getPorts().out.template port<SignalDelayOutputs::z>().connect(&output);
      success *= delay.allocate() == 0 && delay.verify() == 0;
      delay.initialize();
      delay.evaluateResidual();
      success *= output.read() == input && delay.getResidual().getData()[0] == 0;
      delay.tagDifferentiable();
      success *= !delay.tag()[0];
      history.initialize(0, input);
      input = 2.7;
      delay.updateTime(0.1, 0);
      delay.y().getData()[0] = 1.9;
      delay.evaluateResidual();
      success *= std::abs(delay.getResidual().getData()[0] - 0.7072) < 1e-12;
      delay.evaluateJacobian();
      success *= delay.getCooJacobian()->getValues()[0] == -1;
      success *= std::abs(delay.getCooJacobian()->getValues()[1] - 0.768) < 1e-12;

      using AD = DependencyTracking::Variable;
      SignalNode<AD, size_t> ad_source, ad_output;
      AD                     ad_input(2.7, 2);
      ad_source.link(&ad_input, &source_index);
      SignalDelay<AD, size_t> tracked(data, history);
      tracked.getPorts().in.template port<SignalDelayInputs::u>().connect(&ad_source);
      tracked.getPorts().out.template port<SignalDelayOutputs::z>().connect(&ad_output);
      tracked.allocate();
      tracked.initialize();
      tracked.y().getData()[0] = AD(1.9, 0);
      tracked.updateTime(0.1, 0);
      tracked.evaluateResidual();
      const auto& dependencies  = tracked.getResidual().getData()[0].getDependencies();
      success                  *= std::abs(dependencies.at(0) + 1) < 1e-12;
      success                  *= std::abs(dependencies.at(2) - 0.768) < 1e-12;

      history.prepare(0.1, input);
      history.accept();
      delay.evaluateJacobian();
      success      *= delay.getCooJacobian()->getValues()[1] == 0;
      success      *= delay.getCooJacobian()->getColData()[1] == source_index;
      source_index  = INVALID_INDEX<size_t>;
      delay.evaluateJacobian();
      success *= delay.getCooJacobian()->getColData()[1] == 0;
      return success.report(__func__);
    }

    TestOutcome declarations()
    {
      using namespace PhasorDynamics;
      TestStatus                     success = true;
      SignalNodeData<double, size_t> parsed;
      from_json(nlohmann::json::parse(R"({"name":"delayed","signal_id":2,"source":1,"delay":0.01})"), parsed);
      success *= parsed.source == 1 && parsed.delay == 0.01;
      success *= throws([&]
                        { from_json(nlohmann::json::parse(R"({"name":"bad","signal_id":2,"source":1})"), parsed); });
      success *= throws([&]
                        { from_json(nlohmann::json::parse(R"({"name":"bad","signal_id":2,"source":1,"delay":-1})"), parsed); });
      SystemModelData<> data;
      data.signal = {{.name = "root", .signal_id = 0}, {.name = "delayed", .signal_id = 1, .source = 0, .delay = 0.02}, {.name = "chain", .signal_id = 2, .source = 1, .delay = 0.03}, {.name = "alias", .signal_id = 3, .source = 0, .delay = 0}};
      ConstantSignalSourceData<double, size_t> source;
      source.parameters[ConstantSignalSourceParameters::Sr]        = 1.7;
      source.signal_outputs[ConstantSignalSourceSignalOutputs::sr] = 0;
      data.constant_source.push_back(source);
      SystemModel<double, size_t> system(data);
      system.allocate();
      system.initialize();
      success *= system.size() == 2;
      success *= system.getSignalNode(3)->read() == 1.7;
      success *= system.getSignalNode(1)->read() == 1.7;
      success *= throws([&]
                        { system.signalHistory(1); });
      system.addComponent(new InitializationWriter(*system.getSignalNode(1)));
      system.allocate();
      system.initialize();
      success               *= system.getSignalNode(0)->read() == 2.4;
      success               *= system.getSignalNode(1)->read() == 2.4;
      data.signal[0].source  = 2;
      data.signal[0].delay   = 0.01;
      success               *= throws([&]
                        { SystemModel<double, size_t> invalid(data); });
      data.signal[0].source  = 99;
      success               *= throws([&]
                        { SystemModel<double, size_t> invalid(data); });
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
