#include <algorithm>
#include <cmath>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/LogEvaluator.hpp>
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
    class AlgebraicRootEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    protected:
      using NullEvaluator<ScalarT, IdxT>::y_;
      using NullEvaluator<ScalarT, IdxT>::f_;

    public:
      using RealT = typename NullEvaluator<ScalarT, IdxT>::RealT;

      int initialize() override
      {
        NullEvaluator<ScalarT, IdxT>::initialize();
        y_.getData()[0] = 1.0;
        y_.setDataUpdated();
        time_ = 1.0;
        return 0;
      }

      int evaluateResidual() override
      {
        auto*       f = f_.getData();
        const auto* y = y_.getData();
        f[0]          = y[0] * y[0] - static_cast<ScalarT>(time_);
        f_.setDataUpdated();
        return 0;
      }

      void updateTime(RealT t, [[maybe_unused]] RealT a) override
      {
        time_ = t;
      }

    private:
      RealT time_{1.0};
    };

    template <class ScalarT, typename IdxT>
    class DifferentialRampEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    protected:
      using NullEvaluator<ScalarT, IdxT>::y_;
      using NullEvaluator<ScalarT, IdxT>::yp_;
      using NullEvaluator<ScalarT, IdxT>::tag_;
      using NullEvaluator<ScalarT, IdxT>::f_;

    public:
      using RealT = typename NullEvaluator<ScalarT, IdxT>::RealT;

      int initialize() override
      {
        NullEvaluator<ScalarT, IdxT>::initialize();
        y_.getData()[0]  = static_cast<ScalarT>(coordinate_);
        yp_.getData()[0] = 1.0;
        tag_[0]          = true;
        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      int evaluateResidual() override
      {
        auto*       f  = f_.getData();
        const auto* yp = yp_.getData();
        f[0]           = yp[0] - 1.0;
        f_.setDataUpdated();
        return 0;
      }

      void updateTime(RealT coordinate, RealT alpha) override
      {
        coordinate_ = coordinate;
        alpha_      = alpha;
      }

      RealT coordinate() const
      {
        return coordinate_;
      }

      RealT alpha() const
      {
        return alpha_;
      }

    private:
      RealT coordinate_{1.0};
      RealT alpha_{0.0};
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class IdaTests
    {
    public:
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
        TestStatus     success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        ida.setFixedStep(1.0 / n_steps);
        ida.setTolerance(1.0e-6);
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

      TestOutcome algebraicErrorControl()
      {
        Model::AlgebraicRootEvaluator<ScalarT, IdxT> model;

        Ida<ScalarT, IdxT> ida(&model);
        ida.configureSimulation();

        ScalarT max_error = 0.0;
        auto    output_cb = [&](ScalarT t)
        {
          const ScalarT expected = std::sqrt(t);
          max_error              = std::max(max_error, std::abs(model.y().getData()[0] - expected));
        };

        ida.initializeSimulation(1.0, false);
        ida.runSimulation(100.0, (100.0 - 1.0) / 20.0, output_cb);

        TestStatus success  = true;
        success            *= (max_error < 1.0e-3);

        return success.report(__func__);
      }

      TestOutcome logEvaluatorAlgebraic()
      {
        Model::AlgebraicRootEvaluator<ScalarT, IdxT> model;
        Model::LogEvaluator<ScalarT, IdxT>           log_model(model, 1.0);

        log_model.allocate();

        Ida<ScalarT, IdxT> ida(&log_model);
        ida.configureSimulation();

        ScalarT max_error = 0.0;
        auto    output_cb = [&](ScalarT s)
        {
          const ScalarT expected = std::sqrt(std::exp(s));
          max_error              = std::max(max_error, std::abs(log_model.y().getData()[0] - expected));
        };

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(std::log(100.0), std::log(100.0) / 20.0, output_cb);

        TestStatus success  = true;
        success            *= (max_error < 1.0e-3);

        return success.report(__func__);
      }

      TestOutcome logEvaluatorDerivativeScaling()
      {
        Model::DifferentialRampEvaluator<ScalarT, IdxT> model;
        Model::LogEvaluator<ScalarT, IdxT>              log_model(model, 1.0);

        log_model.allocate();

        Ida<ScalarT, IdxT> ida(&log_model);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, false);
        ida.runSimulation(std::log(10.0), std::log(10.0) / 20.0);

        TestStatus success  = true;
        success            *= isEqual(log_model.y().getData()[0], static_cast<ScalarT>(10.0), static_cast<ScalarT>(1.0e-5));

        return success.report(__func__);
      }

      TestOutcome logEvaluatorAlphaScaling()
      {
        Model::DifferentialRampEvaluator<ScalarT, IdxT> model;
        Model::LogEvaluator<ScalarT, IdxT>              log_model(model, 1.0);

        log_model.allocate();
        log_model.initialize();
        log_model.updateTime(std::log(10.0), 30.0);

        TestStatus success  = true;
        success            *= isEqual(model.coordinate(), static_cast<ScalarT>(10.0), static_cast<ScalarT>(1.0e-12));
        success            *= isEqual(model.alpha(), static_cast<ScalarT>(3.0), static_cast<ScalarT>(1.0e-12));

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
