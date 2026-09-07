#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

using AnalysisManager::Sundials::Ida;

namespace GridKit
{
  namespace Model
  {
    template <class ScalarT, typename IdxT>
    class NullEvaluator : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using VectorT = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

      NullEvaluator()
      {
      }

      int allocate() override
      {
        if (!allocated_)
        {
          allocateVectors(size());
          allocated_ = true;
        }
        return 0;
      }

      int initialize()
      {
        if (!allocated_)
        {
          allocate();
        }

        auto* y       = y_.getData();
        auto* yp      = yp_.getData();
        auto* abs_tol = abs_tol_.getData();
        auto* f       = f_.getData();

        y[0]       = 0.0;
        yp[0]      = 0.0;
        tag_       = {false};
        abs_tol[0] = 0.0;
        f[0]       = 0.0;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        abs_tol_.setDataUpdated();
        f_.setDataUpdated();
        return 0;
      }

      IdxT size() override
      {
        return 1;
      }

      IdxT nnz() override
      {
        return 0;
      }

      bool hasJacobian() override
      {
        return false;
      }

      IdxT sizeQuadrature() override
      {
        return 0;
      }

      IdxT sizeParams() override
      {
        return 0;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int setAbsoluteTolerance(RealT rel_tol) override
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      int evaluateResidual() override
      {
        auto*       f = f_.getData();
        const auto* y = y_.getData();
        f[0]          = y[0];
        f_.setDataUpdated();
        return 0;
      }

      int evaluateJacobian() override
      {
        return 0;
      }

      int evaluateIntegrand() override
      {
        return 0;
      }

      int initializeAdjoint() override
      {
        return 0;
      }

      int evaluateAdjointResidual() override
      {
        return 0;
      }

      int evaluateAdjointIntegrand() override
      {
        return 0;
      }

      void updateTime([[maybe_unused]] RealT t, [[maybe_unused]] RealT a) override
      {
      }

      VectorT& y() override
      {
        return y_;
      }

      const VectorT& y() const override
      {
        return y_;
      }

      VectorT& yp() override
      {
        return yp_;
      }

      const VectorT& yp() const override
      {
        return yp_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      VectorT& yB() override
      {
        return yB_;
      }

      const VectorT& yB() const override
      {
        return yB_;
      }

      VectorT& ypB() override
      {
        return ypB_;
      }

      const VectorT& ypB() const override
      {
        return ypB_;
      }

      VectorT& param() override
      {
        return param_;
      }

      const VectorT& param() const override
      {
        return param_;
      }

      VectorT& param_up() override
      {
        return param_up_;
      }

      const VectorT& param_up() const override
      {
        return param_up_;
      }

      VectorT& param_lo() override
      {
        return param_lo_;
      }

      const VectorT& param_lo() const override
      {
        return param_lo_;
      }

      VectorT& getResidual() override
      {
        return f_;
      }

      const VectorT& getResidual() const override
      {
        return f_;
      }

      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* getCsrJacobian() const override
      {
        return csr_jac_;
      }

      VectorT& getIntegrand() override
      {
        return g_;
      }

      const VectorT& getIntegrand() const override
      {
        return g_;
      }

      VectorT& getAdjointResidual() override
      {
        return fB_;
      }

      const VectorT& getAdjointResidual() const override
      {
        return fB_;
      }

      VectorT& getAdjointIntegrand() override
      {
        return gB_;
      }

      const VectorT& getAdjointIntegrand() const override
      {
        return gB_;
      }

      IdxT getIDcomponent()
      {
        return 0;
      }

    protected:
      void allocateVectors(IdxT n)
      {
        y_.resize(n);
        yp_.resize(n);
        f_.resize(n);
        abs_tol_.resize(n);
      }

      VectorT           y_;
      VectorT           yp_;
      std::vector<bool> tag_;
      VectorT           abs_tol_;
      VectorT           f_;
      VectorT           g_;

      VectorT yB_;
      VectorT ypB_;
      VectorT fB_;
      VectorT gB_;

      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* csr_jac_;

      VectorT param_;
      VectorT param_up_;
      VectorT param_lo_;

      bool allocated_{false};
    };

    template <class ScalarT, typename IdxT>
    class AlgebraicErrorControlEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    protected:
      using NullEvaluator<ScalarT, IdxT>::allocated_;
      using NullEvaluator<ScalarT, IdxT>::y_;
      using NullEvaluator<ScalarT, IdxT>::yp_;
      using NullEvaluator<ScalarT, IdxT>::abs_tol_;
      using NullEvaluator<ScalarT, IdxT>::tag_;
      using NullEvaluator<ScalarT, IdxT>::f_;

    public:
      using RealT = typename NullEvaluator<ScalarT, IdxT>::RealT;

      int initialize()
      {
        if (!allocated_)
        {
          this->allocate();
        }

        auto* y       = y_.getData();
        auto* yp      = yp_.getData();
        auto* abs_tol = abs_tol_.getData();
        auto* f       = f_.getData();

        y[0]       = 0.0;
        y[1]       = 0.0;
        yp[0]      = 0.0;
        yp[1]      = 0.0;
        tag_       = {true, false};
        abs_tol[0] = 0.0;
        abs_tol[1] = 0.0;
        f[0]       = 0.0;
        f[1]       = 0.0;
        t_         = 0.0;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        abs_tol_.setDataUpdated();
        f_.setDataUpdated();
        return 0;
      }

      IdxT size() override
      {
        return 2;
      }

      int evaluateResidual() override
      {
        static constexpr RealT OMEGA = 100.0;
        auto*                  f     = f_.getData();
        const auto*            y     = y_.getData();
        const auto*            yp    = yp_.getData();

        f[0] = yp[0];
        f[1] = y[1] - std::sin(OMEGA * t_);
        f_.setDataUpdated();
        return 0;
      }

      void updateTime(RealT t, [[maybe_unused]] RealT a) override
      {
        t_ = t;
      }

    private:
      RealT t_{};
    };

    template <class ScalarT, typename IdxT>
    class ConsistentICTypeEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    protected:
      using NullEvaluator<ScalarT, IdxT>::allocated_;
      using NullEvaluator<ScalarT, IdxT>::y_;
      using NullEvaluator<ScalarT, IdxT>::yp_;
      using NullEvaluator<ScalarT, IdxT>::abs_tol_;
      using NullEvaluator<ScalarT, IdxT>::tag_;
      using NullEvaluator<ScalarT, IdxT>::f_;

    public:
      ConsistentICTypeEvaluator() = default;

      explicit ConsistentICTypeEvaluator(bool steady_state)
        : steady_state_(steady_state)
      {
      }

      int initialize()
      {
        if (!allocated_)
        {
          this->allocate();
        }

        auto* y       = y_.getData();
        auto* yp      = yp_.getData();
        auto* abs_tol = abs_tol_.getData();
        auto* f       = f_.getData();

        y[0]       = 0.0;
        y[1]       = 10.0;
        yp[0]      = steady_state_ ? 0.0 : 2.0; // Purposefully inconsistent guess for non-steady-state case.
        yp[1]      = 0.0;
        tag_       = {true, false};
        abs_tol[0] = 0.0;
        abs_tol[1] = 0.0;
        f[0]       = 0.0;
        f[1]       = 0.0;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        abs_tol_.setDataUpdated();
        f_.setDataUpdated();
        return 0;
      }

      IdxT size() override
      {
        return 2;
      }

      int evaluateResidual() override
      {
        auto*       f  = f_.getData();
        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();

        f[0] = yp[0] + y[0] + y[1] - 1.0;
        f[1] = y[1] - y[0];
        f_.setDataUpdated();
        return 0;
      }

    private:
      bool steady_state_{false};
    };

    // Independent one-state ODE with faults in each evaluator entry point.
    template <class ScalarT, typename IdxT>
    class FaultEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
      using Base = NullEvaluator<ScalarT, IdxT>;

    public:
      using RealT      = typename Base::RealT;
      using VectorT    = typename Base::VectorT;
      using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;

      enum class Operation
      {
        residual,
        jacobian,
        integrand,
        adjoint_residual,
        adjoint_integrand,
        tags,
        tolerance
      };

      struct Failure
      {
        const FaultEvaluator* origin;
        Operation             operation;
      };

      FaultEvaluator()
      {
        jac_.allocateMatrixData(memory::HOST);
        jac_.getRowData()[0] = 0;
        jac_.getRowData()[1] = 1;
        jac_.getColData()[0] = 0;
        jac_.getValues()[0]  = 1.0;
        jac_.setUpdated(memory::HOST);
      }

      int initialize()
      {
        Base::initialize();
        this->y_.setToConst(1.0);
        this->tag_ = {true};
        this->g_.resize(1);
        this->yB_.resize(1);
        this->ypB_.resize(1);
        this->fB_.resize(1);
        this->gB_.resize(1);
        return 0;
      }

      bool hasJacobian() override
      {
        return true;
      }

      IdxT nnz() override
      {
        return 1;
      }

      IdxT sizeQuadrature() override
      {
        return 1;
      }

      IdxT sizeParams() override
      {
        return 1;
      }

      int tagDifferentiable() override
      {
        return fault(Operation::tags);
      }

      int setAbsoluteTolerance(RealT tolerance) override
      {
        const int status = fault(Operation::tolerance);
        return status ? status : Base::setAbsoluteTolerance(tolerance);
      }

      void updateTime(RealT, RealT alpha) override
      {
        alpha_ = alpha;
      }

      int evaluateResidual() override
      {
        if (const int status = fault(Operation::residual))
          return status;
        this->f_.getData()[0] = this->yp_.getData()[0] + this->y_.getData()[0];
        this->f_.setDataUpdated();
        return 0;
      }

      int evaluateJacobian() override
      {
        if (const int status = fault(Operation::jacobian))
          return status;
        jac_.getValues()[0] = 1.0 + alpha_;
        jac_.setUpdated(memory::HOST);
        return 0;
      }

      int evaluateIntegrand() override
      {
        if (const int status = fault(Operation::integrand))
          return status;
        this->g_.getData()[0] = this->y_.getData()[0];
        this->g_.setDataUpdated();
        return 0;
      }

      int initializeAdjoint() override
      {
        this->yB_.setToConst(1.0);
        this->ypB_.setToConst(1.0);
        return 0;
      }

      int evaluateAdjointResidual() override
      {
        if (const int status = fault(Operation::adjoint_residual))
          return status;
        this->fB_.getData()[0] = this->ypB_.getData()[0] - this->yB_.getData()[0];
        this->fB_.setDataUpdated();
        return 0;
      }

      int evaluateAdjointIntegrand() override
      {
        if (const int status = fault(Operation::adjoint_integrand))
          return status;
        this->gB_.getData()[0] = this->yB_.getData()[0];
        this->gB_.setDataUpdated();
        return 0;
      }

      VectorT& getResidual() override
      {
        countRead(Operation::residual);
        return this->f_;
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        countRead(Operation::jacobian);
        return &jac_;
      }

      VectorT& getIntegrand() override
      {
        countRead(Operation::integrand);
        return this->g_;
      }

      VectorT& getAdjointResidual() override
      {
        countRead(Operation::adjoint_residual);
        return this->fB_;
      }

      VectorT& getAdjointIntegrand() override
      {
        countRead(Operation::adjoint_integrand);
        return this->gB_;
      }

      std::optional<Operation> failure;
      int                      status              = 1;
      bool                     throws              = false;
      int                      failures            = 0;
      mutable int              failed_output_reads = 0;

    private:
      int fault(Operation operation)
      {
        if (failure != operation)
          return 0;
        ++failures;
        if (throws)
          throw Failure{this, operation};
        return status;
      }

      void countRead(Operation operation) const
      {
        if (failure == operation)
          ++failed_output_reads;
      }

      mutable CsrMatrixT jac_{1, 1, 1};
      RealT              alpha_{};
    };

    // Method-of-steps oracle: y' = y(t - tau) + A cos(w t), y(t <= 0) = 1.
    template <class ScalarT, typename IdxT>
    class HistoryEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
      using Base = NullEvaluator<ScalarT, IdxT>;

    public:
      using RealT                      = typename Base::RealT;
      static constexpr RealT delay     = 0.2;
      static constexpr RealT amplitude = 10.0;
      static constexpr RealT frequency = 100.0;

      struct Sample
      {
        RealT time, y, yp;
      };

      std::vector<Sample> history;
      std::vector<RealT>  trial_times;
      int                 resets      = 0;
      RealT               reset_value = 0.0;
      RealT               offset      = 0.0;
      RealT               step_limit  = delay;

      int initialize()
      {
        Base::initialize();
        this->y_.setToConst(1.0);
        this->tag_ = {true};
        return 0;
      }

      void resetHistory() override
      {
        reset_value = this->y_.getData()[0];
        history.clear();
        trial_times.clear();
        ++resets;
      }

      void acceptStep(RealT time) override
      {
        history.push_back({time, this->y_.getData()[0], this->yp_.getData()[0]});
      }

      RealT maximumStepSize() const override
      {
        return step_limit;
      }

      void updateTime(RealT time, RealT) override
      {
        time_ = time;
      }

      int evaluateResidual() override
      {
        trial_times.push_back(time_);
        this->f_.getData()[0] = this->yp_.getData()[0] - past(time_ - delay)
                                - amplitude * std::cos(frequency * time_) - offset;
        this->f_.setDataUpdated();
        return 0;
      }

      static RealT exact(RealT time)
      {
        const RealT second = std::max(RealT(0), time - delay);
        return 1.0 + time + amplitude * std::sin(frequency * time) / frequency
               + second * second / 2.0
               + amplitude * (1.0 - std::cos(frequency * second)) / (frequency * frequency);
      }

    private:
      RealT past(RealT time) const
      {
        if (time <= 0.0)
          return 1.0;
        if (history.empty() || time > history.back().time + 1e-12)
          throw std::logic_error("Delay requested unaccepted history");
        const auto next = std::upper_bound(history.begin(), history.end(), time, [](RealT t, const Sample& sample)
                                           { return t < sample.time; });
        if (next == history.end())
          return history.back().y;
        const auto& left = *(next - 1);
        return left.y + (next->y - left.y) * (time - left.time) / (next->time - left.time);
      }

      RealT time_{};
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class IdaTests
    {
    public:
      TestOutcome acceptedHistory()
      {
        using ModelT       = Model::HistoryEvaluator<ScalarT, IdxT>;
        TestStatus success = true;
        ModelT     model;
        model.initialize();
        {
          Ida<ScalarT, IdxT> ida(&model);
          ida.setTolerance(1e-9);
          ida.setMaxSteps(10000);
          ida.configureSimulation();
          success *= model.resets == 1 && model.history.empty();
          ida.initializeSimulation(0.0);
          ida.saveInitialCondition();
          success                      *= model.resets == 2;
          success                      *= model.history.size() == 1;
          success                      *= isEqual(model.history.front().yp, 11.0, 1e-8);
          unsigned            monitors  = 0;
          std::vector<double> monitor_times;
          ida.runSimulation(2.0 * ModelT::delay, 0.013, [&](double time)
                            {
            ++monitors;
            monitor_times.push_back(time);
            success *= isEqual(model.y().getData()[0], ModelT::exact(time), 1e-5); });
          const auto stats  = ida.getStats();
          success          *= model.history.size() == static_cast<size_t>(stats.num_steps_) + 1;
          success          *= stats.num_error_test_fails_ > 0;
          success          *= model.trial_times.size() > model.history.size();
          success          *= monitors == 31;
          success          *= std::any_of(monitor_times.begin(), monitor_times.end(), [&](double t)
                                 { return std::none_of(model.history.begin(), model.history.end(), [&](const auto& sample)
                                                       { return sample.time == t; }); });
          success          *= model.history.back().time == 2.0 * ModelT::delay;
          for (size_t i = 1; i < model.history.size(); ++i)
          {
            const auto& sample  = model.history[i];
            success            *= sample.time > model.history[i - 1].time;
            success            *= sample.time - model.history[i - 1].time <= ModelT::delay + 1e-12;
            success            *= isEqual(sample.y, ModelT::exact(sample.time), 1e-5);
          }

          const auto completed_count  = model.history.size();
          success                    *= ida.runSimulation(2.0 * ModelT::delay) == 0;
          success                    *= model.history.size() == completed_count;

          // A same-time event preserves the left limit and commits a right limit.
          const auto left  = model.history.back();
          const auto count = model.history.size();
          model.offset     = 1.0;
          ida.restartSimulation(left.time);
          success *= model.resets == 2 && model.history.size() == count + 1;
          success *= model.history.back().time == left.time;
          success *= isEqual(model.history.back().y, left.y, 1e-12);
          success *= isEqual(model.history.back().yp, left.yp + 1.0, 1e-6);

          bool rejected_rewind = false;
          try
          {
            ida.restartSimulation(0.0);
          }
          catch (const std::invalid_argument&)
          {
            rejected_rewind = true;
          }
          success *= rejected_rewind && model.history.size() == count + 1;

          // Reuse the configured solver for another study from its saved state.
          model.offset = 0.0;
          ida.getSavedInitialCondition();
          ida.initializeSimulation(0.0);
          success *= model.resets == 3 && model.history.size() == 1;
          success *= model.reset_value == 1.0;
          ida.runSimulation(2.0 * ModelT::delay, 0.013);
          success *= isEqual(model.history.back().y, ModelT::exact(2.0 * ModelT::delay), 1e-5);

          // Loading a new model state has the same explicit fresh-study semantics.
          model.initialize();
          ida.getDefaultInitialCondition();
          ida.initializeSimulation(0.0);
          success *= model.resets == 4 && model.history.size() == 1;
          success *= model.reset_value == 1.0;
        }
        model.initialize();
        model.offset = 0.0;
        {
          Ida<ScalarT, IdxT> ida(&model);
          ida.configureSimulation();
          success *= model.resets == 5 && model.history.empty();
          ida.initializeSimulation(0.0);
          success *= model.resets == 6 && model.history.size() == 1 && model.history.front().time == 0.0;
        }
        return success.report(__func__);
      }

      TestOutcome historyStepLimits()
      {
        TestStatus success = true;
        using RealT        = typename ScalarTraits<ScalarT>::RealT;
        for (const RealT limit : {RealT(0), RealT(-1), std::numeric_limits<RealT>::quiet_NaN()})
        {
          Model::HistoryEvaluator<ScalarT, IdxT> model;
          model.initialize();
          model.step_limit = limit;
          Ida<ScalarT, IdxT> ida(&model);
          bool               caught = false;
          try
          {
            ida.configureSimulation();
          }
          catch (const std::invalid_argument&)
          {
            caught = true;
          }
          success *= caught;
        }
        Model::HistoryEvaluator<ScalarT, IdxT> model;
        model.initialize();
        Ida<ScalarT, IdxT> ida(&model);
        ida.setFixedStep(2.0 * model.delay);
        bool caught = false;
        try
        {
          ida.configureSimulation();
        }
        catch (const std::invalid_argument&)
        {
          caught = true;
        }
        success *= caught;
        // Discrete changes must refresh the bound before a restarted solve.
        Model::HistoryEvaluator<ScalarT, IdxT> restarted_model;
        restarted_model.initialize();
        Ida<ScalarT, IdxT> restarted(&restarted_model);
        restarted.setFixedStep(0.1);
        restarted.configureSimulation();
        restarted.initializeSimulation(0.0);
        restarted_model.step_limit = 0.05;
        caught                     = false;
        try
        {
          restarted.restartSimulation(0.0);
        }
        catch (const std::invalid_argument&)
        {
          caught = true;
        }
        success *= caught;
        return success.report(__func__);
      }

      TestOutcome maximumSteps()
      {
        TestStatus success = true;
        for (const int budget : {-1, 0, 2})
        {
          Model::NullEvaluator<ScalarT, IdxT> model;
          model.initialize();
          Ida<ScalarT, IdxT> ida(&model);
          if (budget >= 0)
            ida.setMaxSteps(static_cast<IdxT>(budget));
          ida.setFixedStep(0.001);
          ida.configureSimulation();
          ida.initializeSimulation(0.0);
          bool caught = false;
          try
          {
            ida.runSimulation(1.0);
          }
          catch (const AnalysisManager::Sundials::SundialsException&)
          {
            caught = true;
          }
          success *= caught;
          success *= ida.getStats().num_steps_ == (budget <= 0 ? 500 : budget);
        }
        return success.report(__func__);
      }

      TestOutcome invalidTimes()
      {
        using RealT                                 = typename ScalarTraits<ScalarT>::RealT;
        TestStatus                          success = true;
        Model::NullEvaluator<ScalarT, IdxT> model;
        model.initialize();
        Ida<ScalarT, IdxT> ida(&model);
        ida.configureSimulation();
        ida.initializeSimulation(0.0);
        const auto rejects = [&](const auto& action)
        {
          bool caught = false;
          try
          {
            action();
          }
          catch (const std::invalid_argument&)
          {
            caught = true;
          }
          success *= caught;
        };
        for (const RealT value : {std::numeric_limits<RealT>::quiet_NaN(), std::numeric_limits<RealT>::infinity()})
        {
          rejects([&]
                  { ida.initializeSimulation(value); });
          rejects([&]
                  { ida.runSimulation(value); });
          rejects([&]
                  { ida.runSimulation(1.0, value); });
        }
        rejects([&]
                { ida.runSimulation(1.0, std::numeric_limits<RealT>::min()); });
        success *= ida.getStats().num_steps_ == 0;
        success *= ida.getInitialTime() == 0.0;
        return success.report(__func__);
      }

      TestOutcome quadratureAndAdjoint()
      {
        TestStatus                           success = true;
        Model::FaultEvaluator<ScalarT, IdxT> model;
        model.initialize();
        Ida<ScalarT, IdxT> ida(&model);
        ida.setTolerance(1e-9);
        ida.setBackwardTolerance(1e-9);
        ida.setQuadratureTolerance(1e-9);
        ida.setBackwardQuadratureTolerance(1e-9);
        ida.configureSimulation();
        ida.configureQuadrature();
        ida.initializeSimulation(0.0);
        ida.initializeQuadrature();
        ida.configureAdjoint();
        ida.initializeAdjoint();
        success *= ida.runForwardSimulation(0.1, 0.013) == 0;
        success *= isEqual(ida.getIntegral()[0], 1.0 - std::exp(-0.1), 1e-7);
        ida.initializeBackwardSimulation(0.1);
        success *= ida.runBackwardSimulation(0.0) == 0;
        success *= isEqual(model.yB().getData()[0], std::exp(-0.1), 1e-7);
        success *= isEqual(ida.getAdjointIntegral()[0], std::exp(-0.1) - 1.0, 1e-7);
        return success.report(__func__);
      }

      TestOutcome evaluationFailures()
      {
        using ModelT       = Model::FaultEvaluator<ScalarT, IdxT>;
        using Operation    = typename ModelT::Operation;
        TestStatus success = true;

        const auto expectFailure = [&](ModelT& model, Operation operation, const auto& action)
        {
          model.failure = operation;
          bool caught   = false;
          try
          {
            action();
          }
          catch (const typename ModelT::Failure& error)
          {
            caught = model.throws && error.origin == &model && error.operation == operation;
          }
          catch (const std::runtime_error& error)
          {
            caught = !model.throws && std::string(error.what()).find("failed with status " + std::to_string(model.status)) != std::string::npos;
          }
          catch (...)
          {
          }
          success *= caught;
          success *= model.failures == 1;
          success *= model.failed_output_reads == 0;
          model.failure.reset();
        };

        // Positive and negative model statuses both mean failure. Exceptions
        // must retain their original type and payload across the C boundary.
        for (const int mode : {1, -1, 0})
        {
          const auto setFault = [mode](ModelT& model)
          {
            model.status = mode;
            model.throws = mode == 0;
            model.initialize();
          };
          for (const auto operation : {Operation::tags, Operation::tolerance})
          {
            ModelT model;
            setFault(model);
            Ida<ScalarT, IdxT> ida(&model);
            expectFailure(model, operation, [&]
                          { ida.configureSimulation(); });
          }

          for (const auto operation : {Operation::residual, Operation::jacobian})
          {
#ifndef GRIDKIT_ENABLE_SUNDIALS_SPARSE
            if (operation == Operation::jacobian)
              continue;
#endif
            for (const bool consistent_ic : {true, false})
            {
              ModelT model;
              setFault(model);
              Ida<ScalarT, IdxT> ida(&model);
              ida.setTolerance(1e-9);
              ida.configureSimulation();
              if (!consistent_ic)
                ida.initializeSimulation(0.0);
              expectFailure(model, operation, [&]
                            {
                if (consistent_ic) ida.initializeSimulation(0.0);
                else ida.runSimulation(0.1); });

              // A handled failure must not poison the next solve.
              model.initialize();
              ida.getDefaultInitialCondition();
              ida.initializeSimulation(0.0);
              ida.runSimulation(0.1);
              success *= isEqual(model.y().getData()[0], std::exp(-0.1), 1e-7);
            }
          }

          for (const bool adjoint : {false, true})
          {
            ModelT model;
            setFault(model);
            Ida<ScalarT, IdxT> ida(&model);
            ida.configureSimulation();
            ida.configureQuadrature();
            ida.initializeSimulation(0.0);
            ida.initializeQuadrature();
            if (adjoint)
            {
              ida.configureAdjoint();
              ida.initializeAdjoint();
            }
            expectFailure(model, Operation::integrand, [&]
                          {
              if (adjoint) ida.runForwardSimulation(0.1);
              else ida.runSimulationQuadrature(0.1); });
          }

          for (const auto operation : {Operation::adjoint_residual, Operation::adjoint_integrand})
          {
            ModelT model;
            setFault(model);
            Ida<ScalarT, IdxT> ida(&model);
            ida.configureSimulation();
            ida.configureQuadrature();
            ida.initializeSimulation(0.0);
            ida.initializeQuadrature();
            ida.configureAdjoint();
            ida.initializeAdjoint();
            ida.runForwardSimulation(0.1);
            ida.initializeBackwardSimulation(0.1);
            expectFailure(model, operation, [&]
                          { ida.runBackwardSimulation(0.0); });
          }
        }
        return success.report(__func__);
      }

      TestOutcome preservesInitialState()
      {
        TestStatus                          success = true;
        Model::NullEvaluator<ScalarT, IdxT> model;
        model.initialize();
        model.y().getData()[0]  = 3.0;
        model.yp().getData()[0] = -2.0;
        Ida<ScalarT, IdxT> ida(&model);
        success *= ida.configureSimulation() == 0;
        success *= model.y().getData()[0] == 3.0;
        success *= model.yp().getData()[0] == -2.0;
        return success.report(__func__);
      }

      TestOutcome callback()
      {
        const unsigned n_steps = 100;
        TestStatus     success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        model.initialize();
        ida.configureSimulation();

        unsigned observed_steps = 0;
        auto     output_cb      = [&]([[maybe_unused]] double t)
        {
          observed_steps++;
        };

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0, 1.0 / n_steps, output_cb);

        success *= (observed_steps == n_steps);

        return success.report(__func__);
      }

      TestOutcome dtMonitorZero()
      {
        TestStatus success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        model.initialize();
        ida.configureSimulation();

        unsigned observed_steps = 0;
        double   observed_t     = 0.0;
        auto     output_cb      = [&](double t)
        {
          observed_steps++;
          observed_t = t;
        };

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0, 0.0, output_cb);

        success *= (observed_steps == 1);
        success *= (observed_t == 1.0);

        return success.report(__func__);
      }

      TestOutcome dtMonitorSuppressesEpsilonFinalStep()
      {
        TestStatus success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        model.initialize();
        ida.configureSimulation();

        unsigned observed_steps = 0;
        double   observed_t     = 0.0;
        auto     output_cb      = [&](double t)
        {
          observed_steps++;
          observed_t = t;
        };

        const double tf = std::nextafter(1.0, 2.0);

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(tf, 0.25, output_cb);

        success *= (observed_steps == 4);
        success *= (observed_t == tf);

        return success.report(__func__);
      }

      TestOutcome fixedStep()
      {
        const unsigned n_steps = 32;
        TestStatus     success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        ida.setFixedStep(1.0 / n_steps);
        ida.setTolerance(1.0e-6);
        model.initialize();
        ida.configureSimulation();

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0);
        auto stats = ida.getStats();

        success *= (stats.num_steps_ == n_steps);

        return success.report(__func__);
      }

      TestOutcome suppressAlgebraicErrors()
      {
        TestStatus success = true;

        const auto countSteps = [](bool suppress_alg)
        {
          Model::AlgebraicErrorControlEvaluator<ScalarT, IdxT> model;

          Ida<ScalarT, IdxT> ida(&model);
          ida.setSuppressAlgebraicErrors(suppress_alg);
          ida.setTolerance(1.0e-6);
          ida.setMaxSteps(10000);
          model.initialize();
          ida.configureSimulation();

          ida.initializeSimulation(0.0, false);
          ida.runSimulation(1.0);

          return ida.getStats().num_steps_;
        };

        const auto unsuppressed_steps = countSteps(false);
        const auto suppressed_steps   = countSteps(true);

        success *= (suppressed_steps < unsuppressed_steps);

        return success.report(__func__);
      }

      TestOutcome consistentICType()
      {
        TestStatus success = true;

        using RealT = typename ScalarTraits<ScalarT>::RealT;

        // If the tolerances are too tight, a finite difference approximation to the Jacobian
        // cannot be generated, and initialization will fail with bad initial guesses.
        static constexpr auto tol = 100.0 * std::numeric_limits<RealT>::epsilon();

        // IdaConsistentICType::YA_YDP with non-steady-state initial guess for the derivatives
        {
          Model::ConsistentICTypeEvaluator<ScalarT, IdxT> model(false);

          Ida<ScalarT, IdxT> ida(&model);
          ida.setConsistentICType(AnalysisManager::Sundials::IdaConsistentICType::YA_YDP);
          ida.setTolerance(tol);
          model.initialize();
          ida.configureSimulation();
          ida.initializeSimulation(0.0);

          success *= isEqual(model.yp().getData()[0], 1.0, tol);
          success *= isEqual(model.yp().getData()[1], 0.0, tol);
          success *= isEqual(model.y().getData()[0], 0.0, tol);
          success *= isEqual(model.y().getData()[1], 0.0, tol);
        }

        // IdaConsistentICType::Y with steady-state initial guess for the derivatives
        {
          Model::ConsistentICTypeEvaluator<ScalarT, IdxT> model(true);

          Ida<ScalarT, IdxT> ida(&model);
          ida.setConsistentICType(AnalysisManager::Sundials::IdaConsistentICType::Y);
          ida.setTolerance(tol);
          model.initialize();
          ida.configureSimulation();
          ida.initializeSimulation(0.0);

          success *= isEqual(model.yp().getData()[0], 0.0, tol);
          success *= isEqual(model.yp().getData()[1], 0.0, tol);
          success *= isEqual(model.y().getData()[0], 0.5, tol);
          success *= isEqual(model.y().getData()[1], 0.5, tol);
        }

        {
          Model::NullEvaluator<ScalarT, IdxT> model;

          Ida<ScalarT, IdxT> ida(&model);
          ida.setConsistentICType(AnalysisManager::Sundials::IdaConsistentICType::Y);
          model.initialize();
          ida.configureSimulation();
          ida.initializeSimulation(0.0);

          success *= isEqual(model.y().getData()[0], 0.0);
        }

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
