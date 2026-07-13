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
      using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using VectorT    = typename Model::Evaluator<ScalarT, IdxT>::VectorT;
      using TagVectorT = typename Model::Evaluator<ScalarT, IdxT>::TagVectorT;

      NullEvaluator()
      {
      }

      int allocate() override
      {
        if (!allocated_)
        {
          allocateVectors(size());
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
        auto* tag     = tag_.getData();
        auto* abs_tol = abs_tol_.getData();
        auto* f       = f_.getData();

        y[0]       = 0.0;
        yp[0]      = 0.0;
        tag[0]     = false;
        abs_tol[0] = 0.0;
        f[0]       = 0.0;
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

      TagVectorT& tag() override
      {
        return tag_;
      }

      const TagVectorT& tag() const override
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

      int setAbsoluteTolerance(RealT rel_tol) override
      {
        std::fill(abs_tol_.getData(), abs_tol_.getData() + abs_tol_.getSize(), rel_tol);
        return 0;
      }

      int evaluateResidual() override
      {
        auto*       f = f_.getData();
        const auto* y = y_.getData();
        f[0]          = y[0];
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
        return nullptr;
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
        y_.allocate(memory::HOST);
        y_.setToZero(memory::HOST);

        yp_.resize(n);
        yp_.allocate(memory::HOST);
        yp_.setToZero(memory::HOST);

        f_.resize(n);
        f_.allocate(memory::HOST);
        f_.setToZero(memory::HOST);

        tag_.resize(n);
        tag_.allocate(memory::HOST);
        tag_.setToConst(false, memory::HOST);

        abs_tol_.resize(n);
        abs_tol_.allocate(memory::HOST);
        abs_tol_.setToZero(memory::HOST);

        allocated_ = true;
      }

      VectorT    y_;
      VectorT    yp_;
      VectorT    f_;
      TagVectorT tag_;
      VectorT    abs_tol_;
      bool       allocated_{false};

      VectorT g_;

      VectorT yB_;
      VectorT ypB_;
      VectorT fB_;
      VectorT gB_;

      VectorT param_;
      VectorT param_up_;
      VectorT param_lo_;
    };

    template <class ScalarT, typename IdxT>
    class AlgebraicErrorControlEvaluator : public NullEvaluator<ScalarT, IdxT>
    {
    public:
      using RealT = typename NullEvaluator<ScalarT, IdxT>::RealT;

      int initialize() override
      {
        if (!this->allocated_)
        {
          this->allocate();
        }

        auto* y       = this->y_.getData();
        auto* yp      = this->yp_.getData();
        auto* tag     = this->tag_.getData();
        auto* abs_tol = this->abs_tol_.getData();
        auto* f       = this->f_.getData();

        y[0]       = 0.0;
        y[1]       = 0.0;
        yp[0]      = 0.0;
        yp[1]      = 0.0;
        tag[0]     = true;
        tag[1]     = false;
        abs_tol[0] = 0.0;
        abs_tol[1] = 0.0;
        f[0]       = 0.0;
        f[1]       = 0.0;
        t_         = 0.0;
        return 0;
      }

      IdxT size() override
      {
        return 2;
      }

      int evaluateResidual() override
      {
        static constexpr RealT OMEGA = 100.0;
        auto*                  f     = this->f_.getData();
        const auto*            y     = this->y_.getData();
        const auto*            yp    = this->yp_.getData();

        f[0] = yp[0];
        f[1] = y[1] - std::sin(OMEGA * t_);
        return 0;
      }

      void updateTime(RealT t, [[maybe_unused]] RealT a) override
      {
        t_ = t;
      }

    private:
      RealT t_{};
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
        ida.runSimulation(1.0, n_steps, output_cb);

        success *= (observed_steps == n_steps);

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
    };
  } // namespace Testing
} // namespace GridKit
