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
    using VectorT = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

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
      return y_.getData()[0];
    }

    const ScalarT& V() const
    {
      return y_.getData()[0];
    }

    // KCL residual accessor
    ScalarT& I()
    {
      return f_.getData()[0];
    }

    const ScalarT& I() const
    {
      return f_.getData()[0];
    }

    // Allocate storage for a single-node voltage and KCL residual
    int allocate()
    {
      size_t size = static_cast<size_t>(size_);

      if (!allocated_)
      {
        allocateVectors(size_);
      }

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
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = V0_;
      yp[0] = 0.0;

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
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    int setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
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
      auto* f = f_.getData();

      f[0] = 0.0;

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

    VectorT           y_;
    VectorT           yp_;
    std::vector<bool> tag_;
    VectorT           abs_tol_;
    VectorT           f_;

    VectorT g_{};
    VectorT param_{};
    VectorT param_up_{};
    VectorT param_lo_{};

    VectorT yB_{};
    VectorT ypB_{};
    VectorT fB_{};
    VectorT gB_{};

    RealT time_{0};
    RealT alpha_{0};

    IdxT max_steps_{0};

    bool allocated_{false};

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

    VectorT& y() final
    {
      return y_;
    }

    const VectorT& y() const final
    {
      return y_;
    }

    VectorT& yp() final
    {
      return yp_;
    }

    const VectorT& yp() const final
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

    VectorT& absoluteTolerance() final
    {
      return abs_tol_;
    }

    const VectorT& absoluteTolerance() const final
    {
      return abs_tol_;
    }

    VectorT& yB() final
    {
      return yB_;
    }

    const VectorT& yB() const final
    {
      return yB_;
    }

    VectorT& ypB() final
    {
      return ypB_;
    }

    const VectorT& ypB() const final
    {
      return ypB_;
    }

    VectorT& param() final
    {
      return param_;
    }

    const VectorT& param() const final
    {
      return param_;
    }

    VectorT& param_up() final
    {
      return param_up_;
    }

    const VectorT& param_up() const final
    {
      return param_up_;
    }

    VectorT& param_lo() final
    {
      return param_lo_;
    }

    const VectorT& param_lo() const final
    {
      return param_lo_;
    }

    VectorT& getResidual() final
    {
      return f_;
    }

    const VectorT& getResidual() const final
    {
      return f_;
    }

    VectorT& getIntegrand() final
    {
      return g_;
    }

    const VectorT& getIntegrand() const final
    {
      return g_;
    }

    VectorT& getAdjointResidual() final
    {
      return fB_;
    }

    const VectorT& getAdjointResidual() const final
    {
      return fB_;
    }

    VectorT& getAdjointIntegrand() final
    {
      return gB_;
    }

    const VectorT& getAdjointIntegrand() const final
    {
      return gB_;
    }

  private:
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

      abs_tol_.resize(n);
      abs_tol_.allocate(memory::HOST);
      abs_tol_.setToZero(memory::HOST);

      allocated_ = true;
    }
  };
} // namespace GridKit
