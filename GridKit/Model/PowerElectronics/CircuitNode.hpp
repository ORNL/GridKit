#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  /**
   * @brief Circuit node representing a connection point.
   */
  template <typename ScalarT, typename IdxT>
  class CircuitNode : public Model::Evaluator<ScalarT, IdxT>
  {
    using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
    using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

  public:
    CircuitNode()
    {
      size_ = 1;
    }

    CircuitNode(ScalarT v0)
      : V0_(v0)
    {
      size_ = 1;
    }

    ~CircuitNode() = default;

    int setNodeID(IdxT id)
    {
      id_ = id;
      return 0;
    }

    IdxT nodeID() const
    {
      return id_;
    }

    // Voltage accessor
    ScalarT& V()
    {
      return y_[0];
    }

    const ScalarT& V() const
    {
      return y_[0];
    }

    // KCL residual accessor
    ScalarT& I()
    {
      return f_[0];
    }

    const ScalarT& I() const
    {
      return f_[0];
    }

    // Allocate storage for a single-node voltage and KCL residual
    int allocate()
    {
      size_t size = static_cast<size_t>(size_);

      y_.resize(size);
      yp_.resize(size);
      f_.resize(size);
      tag_.resize(size);

      variable_indices_[0] = 0;
      residual_indices_[0] = 0;

      return 0;
    }

    /**
     * @brief Initialize node variables
     */
    int initialize()
    {
      y_[0]  = V0_;
      yp_[0] = 0.0;

      return 0;
    }

    /**
     * @brief Node variables are algebraic.
     */
    int tagDifferentiable()
    {
      tag_[0] = false;

      return 0;
    }

    /**
     * @brief Node does not compute residuals, so here we just reset residual values.
     *
     * @warning This implementation assumes node residuals are always evaluated
     * _before_ component model residuals.
     *
     */
    int evaluateResidual()
    {
      f_[0] = 0.0;

      return 0;
    }

    bool hasJacobian() final
    {
      return false;
    }

    /**
     * @brief There is no Jacobian for node variables
     */
    int evaluateJacobian()
    {
      return 0;
    }

    int evaluateIntegrand()
    {
      return 0;
    }

    int initializeAdjoint()
    {
      return 0;
    }

    int evaluateAdjointResidual()
    {
      return 0;
    }

    int evaluateAdjointIntegrand()
    {
      return 0;
    }

  private:
    IdxT    id_{static_cast<IdxT>(-1)};
    IdxT    size_{0};
    IdxT    nnz_{0};
    IdxT    size_quad_{0};
    IdxT    size_opt_{0};
    ScalarT V0_{0.0};

    std::map<IdxT, IdxT> variable_indices_;
    std::map<IdxT, IdxT> residual_indices_;

    std::vector<ScalarT> y_;
    std::vector<ScalarT> yp_;
    std::vector<bool>    tag_;
    std::vector<ScalarT> f_;
    
    MatrixT J_;

    std::vector<ScalarT> g_{};
    std::vector<ScalarT> param_{};
    std::vector<ScalarT> param_up_{};
    std::vector<ScalarT> param_lo_{};

    std::vector<ScalarT> yB_{};
    std::vector<ScalarT> ypB_{};
    std::vector<ScalarT> fB_{};
    std::vector<ScalarT> gB_{};

    RealT time_{0};
    RealT alpha_{0};

    RealT rel_tol_{0};
    RealT abs_tol_{0};

    IdxT max_steps_{0};

  public:
    IdxT size() final
    {
      return size_;
    }

    IdxT nnz() final
    {
      return nnz_;
    }

    IdxT sizeQuadrature() final
    {
      return size_quad_;
    }

    IdxT sizeParams() final
    {
      return size_opt_;
    }

    void updateTime(RealT /* t */, RealT /* a */) final
    {
      // No time to update in node models
    }

    void setTolerances(RealT& rel_tol, RealT& abs_tol) const final
    {
      rel_tol = rel_tol_;
      abs_tol = abs_tol_;
    }

    void setMaxSteps(IdxT& msa) const final
    {
      msa = max_steps_;
    }

    std::vector<ScalarT>& y() final
    {
      return y_;
    }

    const std::vector<ScalarT>& y() const final
    {
      return y_;
    }

    std::vector<ScalarT>& yp() final
    {
      return yp_;
    }

    const std::vector<ScalarT>& yp() const final
    {
      return yp_;
    }

    std::vector<bool>& tag() final
    {
      return tag_;
    }

    const std::vector<bool>& tag() const final
    {
      return tag_;
    }

    std::vector<ScalarT>& yB() final
    {
      return yB_;
    }

    const std::vector<ScalarT>& yB() const final
    {
      return yB_;
    }

    std::vector<ScalarT>& ypB() final
    {
      return ypB_;
    }

    const std::vector<ScalarT>& ypB() const final
    {
      return ypB_;
    }

    std::vector<ScalarT>& param() final
    {
      return param_;
    }

    const std::vector<ScalarT>& param() const final
    {
      return param_;
    }

    std::vector<ScalarT>& param_up() final
    {
      return param_up_;
    }

    const std::vector<ScalarT>& param_up() const final
    {
      return param_up_;
    }

    std::vector<ScalarT>& param_lo() final
    {
      return param_lo_;
    }

    const std::vector<ScalarT>& param_lo() const final
    {
      return param_lo_;
    }

    std::vector<ScalarT>& getResidual() final
    {
      return f_;
    }

    const std::vector<ScalarT>& getResidual() const final
    {
      return f_;
    }

    MatrixT& getJacobian() final
    {
      return J_;
    }

    const MatrixT& getJacobian() const final
    {
      return J_;
    }

    std::vector<ScalarT>& getIntegrand() final
    {
      return g_;
    }

    const std::vector<ScalarT>& getIntegrand() const final
    {
      return g_;
    }

    std::vector<ScalarT>& getAdjointResidual() final
    {
      return fB_;
    }

    const std::vector<ScalarT>& getAdjointResidual() const final
    {
      return fB_;
    }

    std::vector<ScalarT>& getAdjointIntegrand() final
    {
      return gB_;
    }

    const std::vector<ScalarT>& getAdjointIntegrand() const final
    {
      return gB_;
    }
  };
} // namespace GridKit