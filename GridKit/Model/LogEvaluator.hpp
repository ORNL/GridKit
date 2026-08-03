#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Errors.hpp>

namespace GridKit
{
  namespace Model
  {
    template <class ScalarT, typename IdxT>
    class LogEvaluator final : public Evaluator<ScalarT, IdxT>
    {
    public:
      using BaseT      = Evaluator<ScalarT, IdxT>;
      using RealT      = typename BaseT::RealT;
      using VectorT    = typename BaseT::VectorT;
      using CsrMatrixT = typename BaseT::CsrMatrixT;

      LogEvaluator(BaseT& model, RealT omega)
        : model_(model)
      {
        if (omega <= 0.0)
        {
          throw std::invalid_argument("LogEvaluator requires a positive coordinate");
        }

        omega_ = omega;
        s_     = std::log(omega_);
      }

      int allocate() override
      {
        model_.updateTime(omega_, 1.0);
        const int retval = model_.allocate();
        syncSolverDerivativeFromPhysical();
        return retval;
      }

      int initialize() override
      {
        model_.updateTime(omega_, 0.0);
        const int retval = model_.initialize();
        syncSolverDerivativeFromPhysical();
        return retval;
      }

      int tagDifferentiable() override
      {
        return model_.tagDifferentiable();
      }

      int evaluateResidual() override
      {
        return model_.evaluateResidual();
      }

      int evaluateJacobian() override
      {
        return model_.evaluateJacobian();
      }

      int evaluateIntegrand() override
      {
        const int   retval = model_.evaluateIntegrand();
        const auto& g      = model_.getIntegrand();
        const IdxT  n      = g.getSize();
        integrand_log_.resize(n);
        const auto* src = g.getData();
        auto*       dst = integrand_log_.getData();
        for (IdxT i = 0; i < n; ++i)
        {
          dst[i] = src[i] * omega_;
        }
        integrand_log_.setDataUpdated();
        return retval;
      }

      int initializeAdjoint() override
      {
        throw Utilities::NotImplementedError(__func__);
      }

      int evaluateAdjointResidual() override
      {
        throw Utilities::NotImplementedError(__func__);
      }

      int evaluateAdjointIntegrand() override
      {
        throw Utilities::NotImplementedError(__func__);
      }

      IdxT size() override
      {
        return model_.size();
      }

      IdxT nnz() override
      {
        return model_.nnz();
      }

      bool monitoring() const override
      {
        return model_.monitoring();
      }

      void printMonitoredVariables() const override
      {
        model_.printMonitoredVariables();
      }

      const VariableMonitorBase* getMonitor() const override
      {
        return model_.getMonitor();
      }

      void startMonitor() override
      {
        model_.startMonitor();
      }

      void stopMonitor() override
      {
        model_.stopMonitor();
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        return model_.getCsrJacobian();
      }

      bool hasJacobian() override
      {
        return model_.hasJacobian();
      }

      IdxT sizeQuadrature() override
      {
        return model_.sizeQuadrature();
      }

      IdxT sizeParams() override
      {
        return model_.sizeParams();
      }

      void updateTime(RealT s, RealT alpha) override
      {
        s_     = s;
        omega_ = std::exp(s_);

        syncPhysicalDerivativeFromSolver();
        model_.updateTime(omega_, alpha / omega_);
      }

      int setAbsoluteTolerance(RealT rel_tol) override
      {
        return model_.setAbsoluteTolerance(rel_tol);
      }

      VectorT& absoluteTolerance() override
      {
        return model_.absoluteTolerance();
      }

      const VectorT& absoluteTolerance() const override
      {
        return model_.absoluteTolerance();
      }

      VectorT& y() override
      {
        return model_.y();
      }

      const VectorT& y() const override
      {
        return model_.y();
      }

      VectorT& yp() override
      {
        return yp_log_;
      }

      const VectorT& yp() const override
      {
        return yp_log_;
      }

      std::vector<bool>& tag() override
      {
        return model_.tag();
      }

      const std::vector<bool>& tag() const override
      {
        return model_.tag();
      }

      VectorT& yB() override
      {
        return model_.yB();
      }

      const VectorT& yB() const override
      {
        return model_.yB();
      }

      VectorT& ypB() override
      {
        return model_.ypB();
      }

      const VectorT& ypB() const override
      {
        return model_.ypB();
      }

      VectorT& param() override
      {
        return model_.param();
      }

      const VectorT& param() const override
      {
        return model_.param();
      }

      VectorT& param_up() override
      {
        return model_.param_up();
      }

      const VectorT& param_up() const override
      {
        return model_.param_up();
      }

      VectorT& param_lo() override
      {
        return model_.param_lo();
      }

      const VectorT& param_lo() const override
      {
        return model_.param_lo();
      }

      VectorT& getResidual() override
      {
        return model_.getResidual();
      }

      const VectorT& getResidual() const override
      {
        return model_.getResidual();
      }

      VectorT& getIntegrand() override
      {
        return integrand_log_;
      }

      const VectorT& getIntegrand() const override
      {
        return integrand_log_;
      }

      VectorT& getAdjointResidual() override
      {
        return model_.getAdjointResidual();
      }

      const VectorT& getAdjointResidual() const override
      {
        return model_.getAdjointResidual();
      }

      VectorT& getAdjointIntegrand() override
      {
        return model_.getAdjointIntegrand();
      }

      const VectorT& getAdjointIntegrand() const override
      {
        return model_.getAdjointIntegrand();
      }

      RealT omega() const
      {
        return omega_;
      }

      RealT coordinate() const
      {
        return s_;
      }

    private:
      void syncPhysicalDerivativeFromSolver()
      {
        auto&      yp_coordinate = model_.yp();
        const IdxT n             = yp_log_.getSize();
        yp_coordinate.resize(n);

        const auto* src = yp_log_.getData();
        auto*       dst = yp_coordinate.getData();
        for (IdxT i = 0; i < n; ++i)
        {
          dst[i] = src[i] / omega_;
        }
        yp_coordinate.setDataUpdated();
      }

      void syncSolverDerivativeFromPhysical()
      {
        const auto& yp_coordinate = model_.yp();
        const IdxT  n             = yp_coordinate.getSize();
        yp_log_.resize(n);

        const auto* src = yp_coordinate.getData();
        auto*       dst = yp_log_.getData();
        for (IdxT i = 0; i < n; ++i)
        {
          dst[i] = omega_ * src[i];
        }
        yp_log_.setDataUpdated();
      }

      BaseT&  model_;
      RealT   omega_{1.0};
      RealT   s_{0.0};
      VectorT yp_log_;
      VectorT integrand_log_;
    };
  } // namespace Model
} // namespace GridKit
