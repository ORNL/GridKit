#include <cmath>

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

      int initialize() override
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
    class JacobianCountEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    public:
      JacobianCountEvaluator()
      {
        jacobian_.allocateMatrixData(memory::HOST);
        jacobian_.getRowData()[0] = 0;
        jacobian_.getRowData()[1] = 1;
        jacobian_.getColData()[0] = 0;
        jacobian_.getValues()[0]  = 1.0;
        jacobian_.setUpdated(memory::HOST);
        this->csr_jac_ = &jacobian_;
      }

      bool hasJacobian() override
      {
        return true;
      }

      int evaluateJacobian() override
      {
        ++jacobian_calls;
        return 0;
      }

      long int jacobian_calls{0};

    private:
      LinearAlgebra::CsrMatrix<typename NullEvaluator<ScalarT, IdxT>::RealT, IdxT> jacobian_{1, 1, 1};
    };

    /// Nonsingular rotating matrix whose original diagonal pivots vanish at t=0.5.
    template <class ScalarT, typename IdxT>
    class PivotChangeEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
      using RealT = typename NullEvaluator<ScalarT, IdxT>::RealT;

    public:
      PivotChangeEvaluator()
      {
        jacobian_.allocateMatrixData(memory::HOST);
        const IdxT rows[] = {0, 2, 4};
        const IdxT cols[] = {0, 1, 0, 1};
        std::copy_n(rows, 3, jacobian_.getRowData());
        std::copy_n(cols, 4, jacobian_.getColData());
        this->csr_jac_ = &jacobian_;
      }

      IdxT size() override
      {
        return 2;
      }

      IdxT nnz() override
      {
        return 4;
      }

      bool hasJacobian() override
      {
        return true;
      }

      int initialize() override
      {
        this->allocate();
        this->y_.setToZero();
        this->yp_.setToZero();
        this->tag_.assign(2, false);
        return 0;
      }

      void updateTime(RealT t, RealT) override
      {
        time_ = t;
      }

      int evaluateResidual() override
      {
        const ScalarT x = this->y_.getData()[0] - time_ * time_;
        const ScalarT y = this->y_.getData()[1] - time_ * time_ * time_;
        auto*         f = this->f_.getData();
        f[0]            = (1.0 - 2.0 * time_) * x - 2.0 * time_ * y;
        f[1]            = 2.0 * time_ * x + (1.0 - 2.0 * time_) * y;
        this->f_.setDataUpdated();
        return 0;
      }

      int evaluateJacobian() override
      {
        auto* values = jacobian_.getValues();
        values[0] = values[3] = 1.0 - 2.0 * time_;
        values[1]             = -2.0 * time_;
        values[2]             = 2.0 * time_;
        jacobian_.setUpdated(memory::HOST);
        return 0;
      }

    private:
      RealT                                 time_{};
      LinearAlgebra::CsrMatrix<RealT, IdxT> jacobian_{2, 2, 4};
    };

    template <class ScalarT, typename IdxT>
    class MonitoringProbeEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    public:
      explicit MonitoringProbeEvaluator(bool monitoring)
        : monitoring_(monitoring)
      {
      }

      bool monitoring() const override
      {
        ++monitoring_calls_;
        return monitoring_;
      }

      void printMonitoredVariables() const override
      {
        ++print_calls_;
      }

      void resetMonitorCounts()
      {
        monitoring_calls_ = 0;
        print_calls_      = 0;
      }

      std::size_t monitoringCalls() const
      {
        return monitoring_calls_;
      }

      std::size_t printCalls() const
      {
        return print_calls_;
      }

    private:
      bool                monitoring_{};
      mutable std::size_t monitoring_calls_{};
      mutable std::size_t print_calls_{};
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

      int initialize() override
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

      int initialize() override
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
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class IdaTests
    {
    public:
      TestOutcome acceptedStepCallback()
      {
        TestStatus success = true;
        for (const double monitor_interval : {0.0, 0.025})
        {
          std::vector<double>                              baseline;
          std::vector<AnalysisManager::Sundials::IdaStats> baseline_stats;
          for (const bool traced : {false, true})
          {
            Model::ConsistentICTypeEvaluator<ScalarT, IdxT> model;
            Ida<ScalarT, IdxT>                              ida(&model);
            ida.setTolerance(1.0e-8, 1.0e-10);
            ida.configureSimulation();
            ida.initializeSimulation(0.0);
            size_t sample = 0, segment = 0;
            double start = 0.0;
            for (const double end : {0.2, 0.7, 1.0})
            {
              long int                                           accepted = 0;
              double                                             previous = start;
              std::optional<std::function<void(double, double)>> trace;
              if (traced)
              {
                trace = [&](double t, double h)
                {
                  success  *= (std::isfinite(t) && std::isfinite(h) && h > 0 && t > previous);
                  success  *= (std::abs(t - previous - h) < 1.0e-12);
                  previous  = t;
                  ++accepted;
                  success *= (ida.getStats().num_steps_ == accepted);
                };
              }
              auto output = [&](double t)
              {
                const double value  = model.y().getData()[0];
                success            *= (std::abs(value - 0.5 * (1.0 - std::exp(-2.0 * t))) < 1.0e-7);
                if (traced)
                {
                  success *= (sample < baseline.size());
                  if (sample < baseline.size())
                  {
                    success *= (std::abs(value - baseline[sample]) < 1.0e-12);
                  }
                }
                else
                {
                  baseline.push_back(value);
                }
                ++sample;
              };
              ida.runSimulation(end, monitor_interval, output, trace);
              const auto stats = ida.getStats();
              if (traced)
              {
                const auto& expected  = baseline_stats[segment];
                success              *= (accepted == stats.num_steps_ && stats.num_steps_ == expected.num_steps_);
                success              *= (stats.num_residual_evals_ == expected.num_residual_evals_);
                success              *= (stats.num_jacobian_evals_ == expected.num_jacobian_evals_);
                success              *= (stats.num_error_test_fails_ == expected.num_error_test_fails_);
              }
              else
              {
                baseline_stats.push_back(stats);
              }
              ++segment;
              start = end;
              if (end < 1.0)
              {
                ida.initializeSimulation(end);
              }
            }
            success *= (sample == baseline.size());
          }
        }
        return success.report(__func__);
      }

      TestOutcome statisticsAcrossRestarts()
      {
        TestStatus success = true;
#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
        Model::JacobianCountEvaluator<ScalarT, IdxT> model;
        Ida<ScalarT, IdxT>                           ida(&model);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, false);

        AnalysisManager::Sundials::IdaStats totals;
        for (const double end : {0.2, 0.7, 1.0})
        {
          ida.runSimulation(end, 0.0);
          const auto segment  = ida.getStats();
          success            *= (segment.num_steps_ > 0);
          success            *= (segment.num_jacobian_evals_ > 0);
          totals             += segment;
          // Independent oracle: count actual model Jacobian callbacks.
          success            *= (totals.num_jacobian_evals_ == model.jacobian_calls);
          if (end < 1.0)
          {
            ida.initializeSimulation(end, true);
            success *= (ida.getStats().num_steps_ == 0);
          }
        }
        success *= (totals.num_steps_ > ida.getStats().num_steps_);
#endif
        return success.report(__func__);
      }

      TestOutcome callback()
      {
        const unsigned n_steps = 100;
        TestStatus     success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
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

      TestOutcome monitorActivityIsCached()
      {
        TestStatus success = true;

        const auto run = [](bool monitoring)
        {
          Model::MonitoringProbeEvaluator<ScalarT, IdxT> model(monitoring);
          Ida<ScalarT, IdxT>                             ida(&model);
          ida.configureSimulation();
          ida.initializeSimulation(0.0, false);

          model.resetMonitorCounts();
          ida.runSimulation(1.0, 0.25);

          return std::pair(model.monitoringCalls(), model.printCalls());
        };

        const auto [inactive_checks, inactive_prints]  = run(false);
        success                                       *= (inactive_checks == 1);
        success                                       *= (inactive_prints == 0);

        const auto [active_checks, active_prints]  = run(true);
        success                                   *= (active_checks == 1);
        success                                   *= (active_prints == 4);

        return success.report(__func__);
      }

      TestOutcome dtMonitorSuppressesEpsilonFinalStep()
      {
        TestStatus success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
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
        const double   tol     = 1.0e-6;
        TestStatus     success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        ida.setFixedStep(1.0 / n_steps);
        ida.setTolerance(tol);
        ida.configureSimulation();

        // Fixed-step error-test scaling must not affect the tolerance used by
        // the model to construct its absolute-tolerance vector.
        success *= (model.absoluteTolerance().getData()[0] == tol);

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0);
        auto stats = ida.getStats();

        success *= (stats.num_steps_ == n_steps);

        return success.report(__func__);
      }

      TestOutcome changingPivots()
      {
        TestStatus                                 success = true;
        Model::PivotChangeEvaluator<ScalarT, IdxT> model;
        Ida<ScalarT, IdxT>                         ida(&model);
        ida.setFixedStep(0.5);
        ida.setTolerance(1.0e-8);
        ida.configureSimulation();
        try
        {
          ida.initializeSimulation(0.0);
          ida.runSimulation(0.5);
          success *= isEqual(model.y().getData()[0], 0.25, 1.0e-8);
          success *= isEqual(model.y().getData()[1], 0.125, 1.0e-8);
        }
        catch (const AnalysisManager::Sundials::SundialsException&)
        {
          success = false;
        }
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
          ida.configureSimulation();
          ida.initializeSimulation(0.0);

          success *= isEqual(model.y().getData()[0], 0.0);
        }

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
