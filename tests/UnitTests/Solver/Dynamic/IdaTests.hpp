#include "Model/Evaluator.hpp"
#include "Solver/Dynamic/Ida.hpp"
#include "Utilities/TestHelpers.hpp"
#include "Utilities/Testing.hpp"

using AnalysisManager::Sundials::Ida;

namespace GridKit
{
  namespace Model
  {
    template <class ScalarT, typename IdxT>
    class NullEvaluator : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      typedef typename Model::Evaluator<ScalarT, IdxT>::real_type real_type;

      NullEvaluator()
      {
      }

      int allocate() override
      {
        return 0;
      }

      int initialize() override
      {
        y_  = {0};
        yp_ = {0};

        tag_ = {false};

        f_ = {0};
        g_ = {0};
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

      void setTolerances(real_type& rel_tol, real_type& abs_tol) const override
      {
      }

      void setMaxSteps(IdxT& msa) const override
      {
        msa = 2000;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        // f_ = yp_;
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

      void updateTime(real_type t, real_type a) override
      {
      }

      std::vector<ScalarT>& y()
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const
      {
        return y_;
      }

      std::vector<ScalarT>& yp()
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const
      {
        return yp_;
      }

      std::vector<bool>& tag()
      {
        return tag_;
      }

      const std::vector<bool>& tag() const
      {
        return tag_;
      }

      std::vector<ScalarT>& yB()
      {
        return yB_;
      }

      const std::vector<ScalarT>& yB() const
      {
        return yB_;
      }

      std::vector<ScalarT>& ypB()
      {
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const
      {
        return ypB_;
      }

      std::vector<ScalarT>& param()
      {
        return param_;
      }

      const std::vector<ScalarT>& param() const
      {
        return param_;
      }

      std::vector<ScalarT>& param_up()
      {
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const
      {
        return param_up_;
      }

      std::vector<ScalarT>& param_lo()
      {
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const
      {
        return param_lo_;
      }

      std::vector<ScalarT>& getResidual()
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const
      {
        return f_;
      }

      COO_Matrix<ScalarT, IdxT>& getJacobian()
      {
        return jac_;
      }

      const COO_Matrix<ScalarT, IdxT>& getJacobian() const
      {
        return jac_;
      }

      std::vector<ScalarT>& getIntegrand()
      {
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const
      {
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual()
      {
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const
      {
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand()
      {
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const
      {
        return gB_;
      }

      IdxT getIDcomponent()
      {
        return 0;
      }

    protected:
      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;

      std::vector<ScalarT> yB_;
      std::vector<ScalarT> ypB_;
      std::vector<ScalarT> fB_;
      std::vector<ScalarT> gB_;

      COO_Matrix<ScalarT, IdxT> jac_;

      std::vector<ScalarT> param_;
      std::vector<ScalarT> param_up_;
      std::vector<ScalarT> param_lo_;
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class IdaTests
    {
    public:
      TestOutcome test()
      {
        const unsigned n_steps = 100;
        TestStatus     success = true;

        Model::NullEvaluator<ScalarT, IdxT> model;

        Ida<double, size_t> ida(&model);
        ida.configureSimulation();

        unsigned observed_steps = 0;
        auto     output_cb      = [&](double t)
        {
          observed_steps++;
        };

        ida.initializeSimulation(0.0, false);
        ida.runSimulation(1.0, n_steps, output_cb);

        success *= (observed_steps == n_steps);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
