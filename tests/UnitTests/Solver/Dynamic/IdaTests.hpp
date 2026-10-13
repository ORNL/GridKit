#include <array>
#include <cmath>
#include <utility>

#include <idas/idas.h>

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

      TestOutcome nativeStepCheck()
      {
        TestStatus success = true;
#ifndef SUNDIALS_IDAS_STEP_CHECK
        success.skipTest();
#else
        // Exercise the extension directly, independently of the GridKit adapter.
        SUNContext context;
        SUNContext_Create(SUN_COMM_NULL, &context);
        auto y         = N_VNew_Serial(1, context);
        auto yp        = N_VClone(y);
        auto work      = N_VClone(y);
        NV_Ith_S(y, 0) = NV_Ith_S(yp, 0) = 1;
        void* mem                        = IDACreate(context);
        auto  matrix                     = SUNDenseMatrix(1, 1, context);
        auto  linear                     = SUNLinSol_Dense(y, matrix, context);
        IDAInit(mem, [](double, N_Vector y, N_Vector yp, N_Vector r, void*)
                {
                  NV_Ith_S(r, 0) = NV_Ith_S(yp, 0) - NV_Ith_S(y, 0);
                  return 0; },
                0,
                y,
                yp);
        IDASetLinearSolver(mem, linear, matrix);
        IDASStolerances(mem, 1e-8, 1e-10);
        IDASetStopTime(mem, 1);

        struct Check
        {
          void*                                mem;
          N_Vector                             work;
          double                               begin{}, end{}, rejected_begin{};
          int                                  order{}, rejected_order{}, checks{};
          bool                                 rejected{false}, retried{false}, valid{true};
          std::array<std::array<double, 6>, 7> samples{};
        } check{mem, work};

        IDASetStepCheckFn(mem, [](double begin, double end, N_Vector y, N_Vector, double* next, void* data)
                          {
          auto& c = *static_cast<Check*>(data);
          ++c.checks;
          c.begin = begin;
          c.end   = end;
          IDAGetCurrentOrder(c.mem, &c.order);
          if (!c.rejected && begin > 0.1 && c.order >= 3)
          {
            c.rejected       = true;
            c.rejected_begin = begin;
            c.rejected_order = c.order;
            *next            = (end - begin) / 2;
            return 1;
          }
          if (c.rejected && !c.retried)
          {
            c.valid &= begin == c.rejected_begin && c.order == c.rejected_order;
            c.retried = true;
          }
          for (size_t i = 0; i < c.samples.size(); ++i)
          {
            const double t = begin + (end - begin) * static_cast<double>(i) / 6;
            for (int k = 0; k <= c.order; ++k)
            {
              c.valid &= IDAGetDky(c.mem, t, k, c.work) == 0;
              c.samples[i][static_cast<size_t>(k)] = NV_Ith_S(c.work, 0);
            }
          }
          c.valid &= std::abs(c.samples[6][0] - NV_Ith_S(y, 0)) < 1e-12;
          return 0; },
                          &check);
        double t = 0, previous = 0, cap = 0;
        long   steps            = 0;
        int    tightened        = 0;
        double worst_difference = 0;
        while (t < 1 && steps < 1000)
        {
          const int flag  = IDASolve(mem, 1, &t, y, yp, IDA_ONE_STEP);
          success        *= flag >= 0;
          if (flag < 0)
          {
            break;
          }
          ++steps;
          if (cap > 0)
          {
            success *= t - previous <= cap + 1e-12;
          }
          previous = t;
          for (size_t i = 0; i < check.samples.size(); ++i)
          {
            const double q = check.begin + (check.end - check.begin) * static_cast<double>(i) / 6;
            for (int k = 0; k <= check.order; ++k)
            {
              success               *= IDAGetDky(mem, q, k, work) == 0;
              const double expected  = check.samples[i][static_cast<size_t>(k)];
              worst_difference       = std::max(worst_difference, std::abs(NV_Ith_S(work, 0) - expected) / (1 + std::abs(expected)));
              success               *= std::abs(NV_Ith_S(work, 0) - expected) < 1e-9 * (1 + std::abs(expected));
            }
          }
          // Lower a maximum after IDAS has already selected its next step.
          cap = 0;
          if (steps >= 30 && steps < 33)
          {
            double planned;
            IDAGetCurrentStep(mem, &planned);
            cap = std::abs(planned) / 4;
            ++tightened;
          }
          IDASetMaxStep(mem, cap);
        }
        long native_steps;
        IDAGetNumSteps(mem, &native_steps);
        success *= check.valid && check.retried && check.rejected_order >= 3;
        success *= native_steps == steps && check.checks == steps + 1 && tightened == 3;
        success *= std::abs(NV_Ith_S(y, 0) - std::exp(1.0)) < 2e-6;
        std::cout << "native retry retained order " << check.rejected_order
                  << "; relative trial/committed dense difference " << worst_difference << '\n';
        // The native local-error test must also stop at the configured floor.
        IDASetStepCheckFn(mem, nullptr, nullptr);
        NV_Ith_S(y, 0) = NV_Ith_S(yp, 0) = 1;
        IDAReInit(mem, 0, y, yp);
        IDASetMinStep(mem, 0.1);
        IDASStolerances(mem, 1e-12, 1e-14);
        success *= IDASolve(mem, 1, &t, y, yp, IDA_ONE_STEP) == IDA_ERR_FAIL;
        long failures;
        IDAGetNumErrTestFails(mem, &failures);
        success *= failures == 1;
        IDAFree(&mem);
        SUNLinSolFree(linear);
        SUNMatDestroy(matrix);
        N_VDestroy(work);
        N_VDestroy(yp);
        N_VDestroy(y);
        SUNContext_Free(&context);
#endif
        return success.report(__func__);
      }

      TestOutcome controlledStepFailure()
      {
        TestStatus success = true;
        if (!Ida<double, size_t>::supportsStepControl())
        {
          success.skipTest();
          return success.report(__func__);
        }
        for (const auto& [floor, throws] : {std::pair{0.01, false}, std::pair{0.0, false}, std::pair{0.01, true}})
        {
          Model::NullEvaluator<ScalarT, IdxT> model;
          Ida<double, size_t>                 ida(&model);
          unsigned                            accepted = 0, attempts = 0;
          ida.setMinStep(floor);
          ida.setStepControl({[](double, bool) {},
                              [&](double begin, double end, const auto&)
                              {
                                ++attempts;
                                success *= end - begin >= floor;
                                if (throws)
                                {
                                  throw std::runtime_error("step check exception");
                                }
                                return 2.0;
                              },
                              [&]()
                              { ++accepted; },
                              [] {},
                              [] {}});
          ida.configureSimulation();
          ida.initializeSimulation(0, false);
          bool failed = false;
          try
          {
            ida.runSimulation(1);
          }
          catch (const std::exception& error)
          {
            failed = true;
            if (throws)
            {
              success *= std::string(error.what()) == "step check exception";
            }
          }
          success *= failed && accepted == 0 && attempts > 0 && attempts <= 10;
          success *= floor == 0 || attempts == 1;
        }
        return success.report(__func__);
      }

      TestOutcome decreasingStepLimit()
      {
        TestStatus success = true;
        if (!Ida<double, size_t>::supportsStepControl())
        {
          success.skipTest();
          return success.report(__func__);
        }
        Model::NullEvaluator<ScalarT, IdxT> model;
        Ida<double, size_t>                 ida(&model);
        unsigned                            steps = 0;
        double                              limit = 0, largest_step = 0;
        ida.setStepControl({[](double, bool) {},
                            [&](double begin, double end, const auto&)
                            {
                              const double h = end - begin;
                              if (limit > 0)
                              {
                                success *= h <= limit + 1e-12;
                              }
                              largest_step       = std::max(largest_step, h);
                              // Tighten an otherwise unrestricted constant solution,
                              // then release the cap and check that growth resumes.
                              const bool tighten = ++steps >= 4 && steps <= 8;
                              limit              = (tighten ? 0.8 : 2.0) * h;
                              return tighten ? 1.0 : 0.0;
                            },
                            [] {},
                            [] {},
                            [] {}});
        ida.configureSimulation();
        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0);
        success *= steps > 8 && steps < 100;
        success *= largest_step > 0.1;
        success *= ida.getRejectedStepCount() == 0;
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
