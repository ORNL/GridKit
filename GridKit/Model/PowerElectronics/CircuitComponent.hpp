

#pragma once

#include <cassert>
#include <memory>
#include <set>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  /*!
   * @brief Declaration of a CircuitComponent class.
   *
   */
  template <class ScalarT, typename IdxT>
  class CircuitComponent : public Model::Evaluator<ScalarT, IdxT>
  {
  public:
    using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
    using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
    using VectorT    = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

    CircuitComponent() = default;

    /**
     * @note Cannot be marked final, since it is overriden to recurse in the system model.
     */
    void updateTime(RealT t, RealT a) override
    {
      this->time_  = t;
      this->alpha_ = a;
    }

    bool hasJacobian() override
    {
      return true;
    }

    size_t getExternSize()
    {
      return n_extern_;
    }

    size_t getInternalSize()
    {
      return this->n_intern_;
    }

    std::set<size_t> getExternIndices()
    {
      return this->extern_indices_;
    }

    /**
     * @brief Create the mappings from local to global indices
     *
     * @param local_index
     * @param global_index
     * @return int
     */
    int setExternalConnectionNodes(IdxT local_index, IdxT global_index)
    {
      connection_nodes_[static_cast<size_t>(local_index)] = global_index;
      return 0;
    }

    /**
     * @brief Given the location of value in the local vector map to global index
     *
     * f(local_index) = global_index
     *
     * @param local_index index of local value in vector
     * @return IdxT Index of the same value in the global vector
     */
    IdxT getNodeConnection(IdxT local_index) const
    {
      return connection_nodes_[local_index];
    }

    /**
     * @brief Allocates all of the internal buffers for the component.
     * If a components needs a more specialized allocation (such as by having additional internal buffers),
     * it should override this function and then call it in the body to ensure it stays up-to-date with
     * new implementations.
     *
     * @pre \ref nnz_ and \ref size_ must be set. Typically these are set by the child object in its constructor.
     *
     * @return An error code, or 0 if success
     */
    int allocate() override
    {
      jacobian_coo_rows_   = std::make_unique<IdxT[]>(static_cast<size_t>(nnz_));
      jacobian_coo_cols_   = std::make_unique<IdxT[]>(static_cast<size_t>(nnz_));
      jacobian_coo_values_ = std::make_unique<RealT[]>(static_cast<size_t>(nnz_));

      connection_nodes_ = std::make_unique<IdxT[]>(static_cast<size_t>(size_));

      if (!allocated_)
      {
        allocateVectors(size_);
      }

      return 0;
    }

    IdxT* jacobianCooRows()
    {
      return jacobian_coo_rows_.get();
    }

    const IdxT* jacobianCooRows() const
    {
      return jacobian_coo_rows_.get();
    }

    IdxT* jacobianCooCols()
    {
      return jacobian_coo_cols_.get();
    }

    const IdxT* jacobianCooCols() const
    {
      return jacobian_coo_cols_.get();
    }

    RealT* jacobianCooValues()
    {
      return jacobian_coo_values_.get();
    }

    const RealT* jacobianCooValues() const
    {
      return jacobian_coo_values_.get();
    }

    /**
     * @brief Evaluating the residual of a CircuitComponent should be done by evaluating the
     * internal residuals and external residuals. CircuitComponents should overload those
     * functions for their residuals (and the system will call those function instead of this one),
     * so there is no reason to overload this functionality.
     *
     * @return An error code, or 0 is successful.
     */
    int evaluateResidual() final
    {
      if (int err_code = evaluateInternalResidual())
        return err_code;

      return evaluateExternalResidual();
    }

    /**
     * @brief Evaluate all of the residuals of internal variables of the component,
     * modifying \ref f_.
     *
     * @return An error code, or 0 if successful.
     */
    virtual int evaluateInternalResidual() = 0;

    /**
     * @brief Evaluate all of the residuals of external variables of the component,
     * modifying \ref f_
     *
     * @return An error code, or 0 if successful.
     */
    virtual int evaluateExternalResidual() = 0;

    void setInternalPointer(const ScalarT* internals)
    {
      y_int_ = internals;
    }

    void setInternalDerivativePointer(const ScalarT* internals_p)
    {
      yp_int_ = internals_p;
    }

    void setInternalResidualPointer(ScalarT* internal_res)
    {
      f_int_ = internal_res;
    }

  protected:
    /**
     * @brief Reset the Jacobian so it can be constructed. Helper method for \ref setJacValues().
     * Sets \ref current_jac_size_ to 0 so that future calls to `setJacValues()` will override previous values.
     *
     */
    void zeroJacMatrix()
    {
      current_jac_size_ = 0;
    }

    /**
     * @brief Helper method for adding values to the Jacobian. Copies the rows, cols, and vals buffers and appends
     * them to the end of the corresponding Jacobian buffers. Uses \ref current_jac_size_ to tell where the end of
     * the Jacobian currently is.
     *
     * @pre `rows`, `cols`, `vals` must all be the same size
     * @pre \ref allocate() must be called first.
     * @pre Must call \ref zeroJacMatrix() before starting construction of a new Jacobian
     * @pre The must be enough room for the values in the allocated buffers, i.e. `current_jac_size_ + rows.size() <= nnz_`
     */
    void setJacValues(const std::vector<IdxT>& rows, const std::vector<IdxT>& cols, const std::vector<RealT>& vals)
    {
      assert(rows.size() == cols.size());
      assert(rows.size() == vals.size());
      assert(current_jac_size_ + rows.size() <= static_cast<size_t>(nnz_));

      for (size_t i = 0; i < rows.size(); i++)
      {
        jacobian_coo_rows_[current_jac_size_]   = rows[i];
        jacobian_coo_cols_[current_jac_size_]   = cols[i];
        jacobian_coo_values_[current_jac_size_] = vals[i];

        current_jac_size_++;
      }
    }

  public:
    IdxT size() final
    {
      return size_;
    }

    IdxT size() const
    {
      return size_;
    }

    IdxT nnz() final
    {
      return nnz_;
    }

    IdxT nnz() const
    {
      return nnz_;
    }

    IdxT sizeQuadrature() final
    {
      return size_quad_;
    }

    IdxT sizeQuadrature() const
    {
      return size_quad_;
    }

    IdxT sizeParams() final
    {
      return size_opt_;
    }

    IdxT sizeParams() const
    {
      return size_opt_;
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

    //@todo Fix ID naming
    IdxT getIDcomponent() const
    {
      return idc_;
    }

  protected:
    /**
     * @brief Allocate state and residual storage owned by this component.
     */
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

    size_t                  n_extern_;
    size_t                  n_intern_;
    std::set<IdxT>          extern_indices_;
    ///@todo may want to replace the mapping of connection_nodes to Node objects instead of IdxT. Allows for container free setup
    std::unique_ptr<IdxT[]> connection_nodes_;

  protected:
    IdxT size_{0};
    IdxT nnz_{0};
    IdxT size_quad_{0};
    IdxT size_opt_{0};

    // COO Jacobian buffers
    std::unique_ptr<IdxT[]>  jacobian_coo_rows_;
    std::unique_ptr<IdxT[]>  jacobian_coo_cols_;
    std::unique_ptr<RealT[]> jacobian_coo_values_;

    /// The number of non-zero elements currently inserted into the Jacobian. See \ref setJacValues()
    size_t current_jac_size_{0};

    /// @brief A pointer to the internal variables of this component.
    const ScalarT* y_int_;
    /// @brief A pointer to the internal derivatives of this component.
    const ScalarT* yp_int_;
    /// @brief A pointer to the internal residuals of this component
    ScalarT*       f_int_;

    VectorT           y_;
    VectorT           yp_;
    std::vector<bool> tag_;
    VectorT           abs_tol_;
    VectorT           f_;

    VectorT g_;

    VectorT yB_;
    VectorT ypB_;
    VectorT fB_;
    VectorT gB_;

    VectorT param_;
    VectorT param_up_;
    VectorT param_lo_;

    RealT time_;
    RealT alpha_;

    IdxT max_steps_;

    IdxT idc_;

    bool allocated_{false};
  };

} // namespace GridKit
