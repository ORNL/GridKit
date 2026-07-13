#pragma once

#include <functional>

#include "Evaluator.hpp"

namespace GridKit
{
  namespace Model
  {

    template <class ScalarT, typename IdxT>
    class PartitionEvaluator : public Evaluator<ScalarT, IdxT>
    {
    public:
      using Base       = Evaluator<ScalarT, IdxT>;
      using RealT      = typename Base::RealT;
      using MatrixT    = typename Base::MatrixT;
      using CsrMatrixT = typename Base::CsrMatrixT;

      using ForcingFunc = std::function<ScalarT(RealT)>;

    private:
      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<ScalarT> yB_;
      std::vector<ScalarT> ypB_;
      std::vector<ScalarT> param_;
      std::vector<ScalarT> param_up_;
      std::vector<ScalarT> param_lo_;
      std::vector<ScalarT> residual_;
      std::vector<ScalarT> integrand_;
      std::vector<ScalarT> adj_integrand_;
      std::vector<ScalarT> adj_residual_;
      std::vector<bool>    tag_;

      MatrixT jacobian_;

      std::vector<ScalarT>     coupling_;
      std::vector<ForcingFunc> forcing_;

    public:
      int allocate() override
      {
        return 0;
      }

      int initialize() override
      {
        return 0;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
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

      IdxT size() override
      {
        return 0;
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

      void updateTime(RealT, RealT) override
      {
      }

      void setTolerances(RealT&, RealT&) const override
      {
      }

      void setMaxSteps(IdxT&) const override
      {
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
        return residual_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return residual_;
      }

      MatrixT& getJacobian() override
      {
        return jacobian_;
      }

      const MatrixT& getJacobian() const override
      {
        return jacobian_;
      }

      std::vector<ScalarT>& getIntegrand() override
      {
        return integrand_;
      }

      const std::vector<ScalarT>& getIntegrand() const override
      {
        return integrand_;
      }

      std::vector<ScalarT>& getAdjointResidual() override
      {
        return adj_residual_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override
      {
        return adj_residual_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override
      {
        return adj_integrand_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        return adj_integrand_;
      }

      void addForcing(const ForcingFunc& f)
      {
        forcing_.push_back(f);
      }

      IdxT getCouplingSize() const
      {
        return static_cast<IdxT>(coupling_.size());
      }

      IdxT getStateSize() const
      {
        return static_cast<IdxT>(y_.size());
      }

      std::vector<ScalarT>& getCoupling()
      {
        return coupling_;
      }

      const std::vector<ScalarT>& getCoupling() const
      {
        return coupling_;
      }
    };

  } // namespace Model
} // namespace GridKit
