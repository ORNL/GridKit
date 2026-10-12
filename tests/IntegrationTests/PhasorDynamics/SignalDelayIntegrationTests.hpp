#pragma once

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayControl.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  class SignalDelayIntegrationTests
  {
    using System = PhasorDynamics::SystemModel<double, size_t>;
    using Solver = AnalysisManager::Sundials::Ida<double, size_t>;
    using Node   = PhasorDynamics::SignalNode<double, size_t>;

    // x' = -gamma*x + k*x(t-tau), with exact smooth solution exp(lambda*t).
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

    struct Result
    {
      double error{0};
      double maximum_step{0};
      double final_value{0};
      size_t outputs{0};
      size_t rejected{0};
      bool   rollback_ok{true};
    };

    Result solve(double tolerance, bool compress, bool force_rejection = false, bool restart = false)
    {
      using namespace PhasorDynamics;
      constexpr double  tau = 0.001;
      SystemModelData<> data;
      data.signal = {{.name = "x", .signal_id = 0}, {.name = "short", .signal_id = 1, .source = 0, .delay = tau}, {.name = "long", .signal_id = 2, .source = 1, .delay = 0.199}};
      System system(data);
      system.addComponent(new Feedback(*system.getSignalNode(0), *system.getSignalNode(1), tau));
      system.signalHistory(0).setPrehistory([](double t)
                                            { return std::exp(-0.4 * t); });
      system.allocate();
      system.stopMonitor();
      system.setSignalHistoryTolerance(tolerance, tolerance * 0.01, compress);
      Solver solver(&system);
      solver.setTolerance(tolerance, tolerance * 0.01);
      solver.setMaxSteps(100000);
      Result result;
      bool   rejected        = false;
      double candidate_begin = 0, candidate_end = 0;
      size_t output_count_before_rejection = 0;
      solver.setStepControl({[&](double t, bool consistent)
                             { system.initializeSignalHistory(t, consistent); },
                             [&](double begin, double end, const Solver::DenseOutput& dense)
                             {
                               candidate_begin    = begin;
                               candidate_end      = end;
                               const double error = system.assessSignalHistory(begin, end, dense);
                               if (force_rejection && !rejected && end > 0.1)
                               {
                                 rejected                      = true;
                                 output_count_before_rejection = result.outputs;
                                 return 2.0;
                               }
                               return error;
                             },
                             [&]()
                             {
                               system.acceptSignalHistory();
                               result.maximum_step = std::max(result.maximum_step, candidate_end - candidate_begin);
                             },
                             [&]()
                             { system.finishSignalHistory(); },
                             [&]()
                             {
                               system.discardSignalHistory();
                               result.rollback_ok &= system.signalHistory(0).time() == candidate_begin;
                               if (force_rejection && rejected)
                               {
                                 result.rollback_ok &= result.outputs == output_count_before_rejection;
                               }
                             }});
      solver.configureSimulation();
      solver.initializeSimulation(0);
      auto observe = [&](double t)
      {
        ++result.outputs;
        result.error = std::max(result.error, std::abs(system.getSignalNode(0)->read() - std::exp(-0.4 * t)));
        result.error = std::max(result.error, std::abs(system.getSignalNode(1)->read() - std::exp(-0.4 * (t - tau))));
        result.error = std::max(result.error, std::abs(system.getSignalNode(2)->read() - std::exp(-0.4 * (t - 0.2))));
      };
      if (restart)
      {
        solver.runSimulation(0.5, 0.025, observe);
        const auto count = system.signalHistory(0).segmentCount();
        solver.initializeSimulation(0.5);
        result.rollback_ok &= system.signalHistory(0).segmentCount() == count;
      }
      solver.runSimulation(1, 0.025, observe);
      result.rejected    = solver.getRejectedStepCount();
      result.final_value = system.getSignalNode(0)->read();
      std::cout << "delay tolerance=" << tolerance << " error=" << result.error
                << " max h/tau=" << result.maximum_step / tau << " rejected=" << result.rejected << '\n';
      return result;
    }

  public:
    TestOutcome prescribedSource()
    {
      using namespace PhasorDynamics;
      TestStatus        success = true;
      SystemModelData<> data;
      data.signal = {{.name = "u", .signal_id = 0},
                     {.name = "z", .signal_id = 1, .source = 0, .delay = 2.0}};
      FunctionSignalSourceData<double, size_t> source;
      source.parameters[FunctionSignalSourceParameters::Fr]        = std::string("cos");
      source.parameters[FunctionSignalSourceParameters::Fi]        = std::string("sin");
      source.signal_outputs[FunctionSignalSourceSignalOutputs::sr] = 0;
      data.function_source.push_back(source);
      System system(data);
      system.allocate();
      system.stopMonitor();
      Solver solver(&system);
      solver.setTolerance(1e-4, 1e-6);
      solver.setMaxSteps(100000);
      configureSignalDelays(system, solver, 1e-10, 1e-12, false);
      solver.configureSimulation();
      solver.initializeSimulation(0);
      solver.runSimulation(1);
      // The output is still in constant prehistory. Its source must nevertheless
      // have accurate stored history for future delayed reads.
      success *= std::abs(system.getSignalNode(1)->read() - 1) < 1e-12;
      for (int i = 1; i < 100; ++i)
      {
        const double q  = i * 0.01;
        success        *= std::abs(system.signalHistory(0).sample(1, 1 - q).constant - std::cos(q)) < 2e-9;
      }
      success *= solver.getRejectedStepCount() > 0;
      // A fresh run may start at a nonzero time with explicitly supplied past.
      solver.getDefaultInitialCondition();
      system.signalHistory(0).setPrehistory([](double t)
                                            { return std::cos(t); });
      solver.initializeSimulation(2);
      success *= system.signalHistory(0).time() == 2;
      success *= system.signalHistory(0).segmentCount() == 0;
      solver.runSimulation(2.1);
      success *= std::abs(system.getSignalNode(1)->read() - std::cos(0.1)) < 1e-5;
      return success.report(__func__);
    }

    TestOutcome feedbackAndConvergence()
    {
      TestStatus success  = true;
      const auto coarse   = solve(1e-4, false);
      const auto fine     = solve(1e-7, false);
      const auto compact  = solve(1e-7, true);
      success            *= fine.error < coarse.error;
      success            *= fine.error < 5e-6;
      success            *= compact.error < 5e-6;
      success            *= coarse.maximum_step > 0.001 && fine.maximum_step > 0.001;
      success            *= std::abs(fine.final_value - compact.final_value) < 2e-6;
      success            *= fine.outputs == 40;
      return success.report(__func__);
    }

    TestOutcome rollbackAndContinuation()
    {
      TestStatus success    = true;
      const auto rejected   = solve(1e-6, false, true);
      const auto continued  = solve(1e-6, false, false, true);
      success              *= rejected.rejected > 0 && rejected.rollback_ok;
      success              *= rejected.outputs == 40 && rejected.error < 5e-5;
      success              *= continued.outputs == 40 && continued.error < 5e-5 && continued.rollback_ok;
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
