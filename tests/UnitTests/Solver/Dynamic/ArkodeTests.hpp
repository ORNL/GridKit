#include <cmath>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/ArkStep.hpp>
#include <GridKit/Solver/Dynamic/ErkStep.hpp>
#include <GridKit/Solver/Dynamic/SundialsException.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

using AnalysisManager::Sundials::ArkStep;
using AnalysisManager::Sundials::ErkStep;
using AnalysisManager::Sundials::SundialsException;

namespace GridKit
{
  namespace Model
  {
    template <class ScalarT, typename IdxT>
    class LinearOdeEvaluator : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using Base       = Model::Evaluator<ScalarT, IdxT>;
      using RealT      = typename Base::RealT;
      using MatrixT    = typename Base::MatrixT;
      using CsrMatrixT = typename Base::CsrMatrixT;

      LinearOdeEvaluator(ScalarT mass, ScalarT lambda)
        : mass_(mass),
          lambda_(lambda),
          csr_jac_(1, 1, 1)
      {
        IdxT  rows[2] = {0, 1};
        IdxT  cols[1] = {0};
        RealT vals[1] = {static_cast<RealT>(lambda_)};
        csr_jac_.copyDataFrom(rows, cols, vals, LinearAlgebra::memory::HOST, LinearAlgebra::memory::HOST);
      }

      int allocate() override
      {
        return 0;
      }

      int initialize() override
      {
        y_       = {1.0};
        yp_      = {0.0};
        tag_     = {true};
        abs_tol_ = {0.0};
        f_       = {0.0};
        g_.clear();
        yB_.clear();
        ypB_.clear();
        fB_.clear();
        gB_.clear();
        param_.clear();
        param_up_.clear();
        param_lo_.clear();
        return 0;
      }

      int tagDifferentiable() override
      {
        tag_ = {true};
        return 0;
      }

      int setAbsoluteTolerance(RealT rel_tol) override
      {
        abs_tol_ = {rel_tol};
        return 0;
      }

      int evaluateResidual() override
      {
        f_[0] = lambda_ * y_[0] - mass_ * yp_[0];
        return 0;
      }

      int evaluateJacobian() override
      {
        RealT vals[1] = {static_cast<RealT>(lambda_) - alpha_ * static_cast<RealT>(mass_)};
        csr_jac_.copyValues(vals, LinearAlgebra::memory::HOST, LinearAlgebra::memory::HOST);
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

      IdxT size() override
      {
        return 1;
      }

      IdxT nnz() override
      {
        return 1;
      }

      bool hasJacobian() override
      {
        return true;
      }

      IdxT sizeQuadrature() override
      {
        return 0;
      }

      IdxT sizeParams() override
      {
        return 0;
      }

      void updateTime(RealT, RealT a) override
      {
        alpha_ = a;
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        return &csr_jac_;
      }

      std::vector<ScalarT>& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      std::vector<ScalarT>& y() override
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const override
      {
        return y_;
      }

      std::vector<ScalarT>& yp() override
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const override
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

      std::vector<ScalarT>& yB() override
      {
        return yB_;
      }

      const std::vector<ScalarT>& yB() const override
      {
        return yB_;
      }

      std::vector<ScalarT>& ypB() override
      {
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const override
      {
        return ypB_;
      }

      std::vector<ScalarT>& param() override
      {
        return param_;
      }

      const std::vector<ScalarT>& param() const override
      {
        return param_;
      }

      std::vector<ScalarT>& param_up() override
      {
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const override
      {
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() override
      {
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const override
      {
        return param_lo_;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      MatrixT& getJacobian() override
      {
        return jac_;
      }

      const MatrixT& getJacobian() const override
      {
        return jac_;
      }

      std::vector<ScalarT>& getIntegrand() override
      {
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const override
      {
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() override
      {
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override
      {
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override
      {
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        return gB_;
      }

    private:
      ScalarT mass_;
      ScalarT lambda_;
      RealT   alpha_{};

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> abs_tol_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;
      std::vector<ScalarT> yB_;
      std::vector<ScalarT> ypB_;
      std::vector<ScalarT> fB_;
      std::vector<ScalarT> gB_;
      std::vector<ScalarT> param_;
      std::vector<ScalarT> param_up_;
      std::vector<ScalarT> param_lo_;

      MatrixT            jac_;
      mutable CsrMatrixT csr_jac_;
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ArkodeTests
    {
      using RealT = typename Model::Evaluator<ScalarT, IdxT>::RealT;

    public:
      TestOutcome arkStepIdentityMass()
      {
        return runArkStep(1.0, __func__);
      }

      TestOutcome arkStepNonIdentityMass()
      {
        return runArkStep(2.0, __func__);
      }

      TestOutcome erkStepForwardEuler()
      {
        const unsigned n_steps = 64;
        const RealT    h       = 1.0 / static_cast<RealT>(n_steps);
        const RealT    lambda  = -0.5;
        TestStatus     success = true;

        Model::LinearOdeEvaluator<ScalarT, IdxT> model(1.0, lambda);
        ErkStep<ScalarT, IdxT>                   solver(&model);
        solver.setOrder(1);
        solver.setFixedStep(h);
        solver.setTolerance(1.0e-12);
        solver.configureSimulation();
        solver.initializeSimulation(0.0);
        solver.runSimulation(1.0);

        const RealT expected  = std::pow(1.0 + h * lambda, static_cast<RealT>(n_steps));
        success              *= isEqual(model.y()[0], expected, 1.0e-10);

        return success.report(__func__);
      }

      TestOutcome erkStepRejectsNonIdentityMass()
      {
        TestStatus success = true;
        bool       threw   = false;

        try
        {
          Model::LinearOdeEvaluator<ScalarT, IdxT> model(2.0, -0.5);
          ErkStep<ScalarT, IdxT>                   solver(&model);
          solver.configureSimulation();
        }
        catch (const SundialsException&)
        {
          threw = true;
        }

        success *= threw;
        return success.report(__func__);
      }

      TestOutcome callback()
      {
        const int  nout    = 12;
        TestStatus success = true;

        Model::LinearOdeEvaluator<ScalarT, IdxT> model(1.0, -0.5);
        ArkStep<ScalarT, IdxT>                   solver(&model);
        solver.configureSimulation();

        int  observed = 0;
        auto step_cb  = [&]([[maybe_unused]] RealT t)
        {
          ++observed;
        };

        solver.initializeSimulation(0.0);
        solver.runSimulation(1.0, nout, step_cb);

        success *= (observed == nout);
        return success.report(__func__);
      }

    private:
      TestOutcome runArkStep(ScalarT mass, const char* name)
      {
        const RealT lambda  = -0.5;
        TestStatus  success = true;

        Model::LinearOdeEvaluator<ScalarT, IdxT> model(mass, lambda);
        ArkStep<ScalarT, IdxT>                   solver(&model);
        solver.setTolerance(1.0e-9);
        solver.configureSimulation();
        solver.initializeSimulation(0.0);
        solver.runSimulation(1.0, 8);

        const RealT expected  = std::exp(lambda / static_cast<RealT>(mass));
        success              *= isEqual(model.y()[0], expected, 1.0e-6);

        return success.report(name);
      }
    };
  } // namespace Testing
} // namespace GridKit
