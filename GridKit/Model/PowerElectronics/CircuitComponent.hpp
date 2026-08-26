

#pragma once

#include <cassert>
#include <memory>
#include <set>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/ExternalConnection.hpp>

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

    CircuitComponent(const CircuitComponent& other)
      : n_extern_(other.n_extern_),
        n_intern_(other.n_intern_),
        extern_indices_(other.extern_indices_),
        size_(other.size_),
        nnz_(other.nnz_),
        size_quad_(other.size_quad_),
        size_opt_(other.size_opt_),
        current_jac_size_(other.current_jac_size_),

        // These pointers refer to storage supplied by a parent system.
        // The copied component must be connected to its own storage later.
        y_int_(nullptr),
        yp_int_(nullptr),
        f_int_(nullptr),

        tag_(other.tag_),
        time_(other.time_),
        alpha_(other.alpha_),
        max_steps_(other.max_steps_),
        idc_(other.idc_),
        allocated_(other.allocated_)
    {
      /*
       * VectorT disables its normal copy constructor and copy-assignment
       * operator. Use its provided copyFromExternal() operation to perform
       * an independent copy of the vector data.
       */
      auto copyVector = [](VectorT& destination, const VectorT& source)
      {
        const IdxT source_size = source.getSize();

        if (source_size == 0)
        {
          return;
        }

        destination.resize(source_size);
        destination.copyFromExternal(source);
      };

      /*
       * Deep-copy the local-to-global connection mapping.
       */
      if (other.connection_nodes_)
      {
        connection_nodes_ = std::make_unique<IdxT[]>(static_cast<size_t>(size_));

        for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
        {
          connection_nodes_[i] = other.connection_nodes_[i];
        }
      }

      /*
       * Deep-copy the COO Jacobian row indices.
       */
      if (other.jacobian_coo_rows_)
      {
        jacobian_coo_rows_ = std::make_unique<IdxT[]>(static_cast<size_t>(nnz_));

        for (size_t i = 0; i < static_cast<size_t>(nnz_); ++i)
        {
          jacobian_coo_rows_[i] = other.jacobian_coo_rows_[i];
        }
      }

      /*
       * Deep-copy the COO Jacobian column indices.
       */
      if (other.jacobian_coo_cols_)
      {
        jacobian_coo_cols_ = std::make_unique<IdxT[]>(static_cast<size_t>(nnz_));

        for (size_t i = 0; i < static_cast<size_t>(nnz_); ++i)
        {
          jacobian_coo_cols_[i] = other.jacobian_coo_cols_[i];
        }
      }

      /*
       * Deep-copy the COO Jacobian values.
       */
      if (other.jacobian_coo_values_)
      {
        jacobian_coo_values_ = std::make_unique<RealT[]>(static_cast<size_t>(nnz_));

        for (size_t i = 0; i < static_cast<size_t>(nnz_); ++i)
        {
          jacobian_coo_values_[i] = other.jacobian_coo_values_[i];
        }
      }

      if (size_ > 0)
      {
        y_ext_  = std::make_unique<const ScalarT*[]>(static_cast<size_t>(size_));
        yp_ext_ = std::make_unique<const ScalarT*[]>(static_cast<size_t>(size_));
        f_ext_  = std::make_unique<ScalarT*[]>(static_cast<size_t>(size_));

        for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
        {
          y_ext_[i]  = nullptr;
          yp_ext_[i] = nullptr;
          f_ext_[i]  = nullptr;
        }
      }

      // State, state derivative, residual, and absolute tolerance.
      copyVector(y_, other.y_);
      copyVector(yp_, other.yp_);
      copyVector(f_, other.f_);
      copyVector(abs_tol_, other.abs_tol_);
      copyVector(g_, other.g_);
      copyVector(yB_, other.yB_);
      copyVector(ypB_, other.ypB_);
      copyVector(fB_, other.fB_);
      copyVector(gB_, other.gB_);
      copyVector(param_, other.param_);
      copyVector(param_up_, other.param_up_);
      copyVector(param_lo_, other.param_lo_);
    }

    virtual CircuitComponent<ScalarT, IdxT>* clone() const
    {
      return nullptr;
    }

    virtual bool isCloneable() const
    {
      return false;
    }

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

    std::set<IdxT> getExternIndices()
    {
      return this->extern_indices_;
    }

    /**
     * @brief Create the mappings from local to global indices for an internal variable.
     * Used for constructing system Jacobians \see connection_nodes_.
     *
     * @param local_index The index of the local variable
     * @param global_index The index of the corresponding system variable.
     *
     * @pre `local_index` *must* be the index of an internal variable. Using this method for
     * an external variable will not properly setup the data pointers for that variable.
     */
    int setInternalConnectionNodes(size_t local_index, IdxT global_index)
    {
      assert(!extern_indices_.contains(static_cast<IdxT>(local_index)));
      setConnectionNodes(local_index, global_index);
      return 0;
    }

    /**
     * @brief Create the mappings from local to global indices for an external variable.
     * External variables need extra information than internal variables - their data
     * pointers \ref y_ext_, \ref yp_ext_, and \ref f_ext_.
     *
     * @param local_index The index of the local variable
     * @param connection The necessary connection information for the variable
     *
     * @pre `local_index` *must* be the index of an external variable. As of now, using this method
     * to set information for a local variable will silently discard the unnecessary information, but
     * this may change in the future.
     */
    int setExternalConnectionNodes(size_t local_index, ExternalConnection<ScalarT, IdxT> connection)
    {
      assert(extern_indices_.contains(local_index));
      y_ext_[local_index]  = connection.y_;
      yp_ext_[local_index] = connection.yp_;
      f_ext_[local_index]  = connection.f_;
      setConnectionNodes(local_index, connection.idx_);
      return 0;
    }

    /**
     * @brief Update the connection index for a variable.
     *
     * Sets only the connection index without modifying the variable's
     * internal/external classification or its associated data pointers.
     *
     * @param local_index Index of the local variable.
     * @param connection_index New connection index for the variable.
     *
     * @return int 0 if successful.
     */
    int setConnectionNodes(size_t local_index, IdxT connection_index)
    {
      connection_nodes_[local_index] = connection_index;
      return 0;
    }

    /**
     * @brief Given the location of value in the local vector map to global index
     *
     * f(local_index) = global_index
     *
     * @param local_index index of local value in vector
     * @return size_t Index of the same value in the global vector
     */
    IdxT getNodeConnection(size_t local_index) const
    {
      return connection_nodes_[local_index];
    }

    int initialize() override
    {
      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
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

      y_ext_  = std::make_unique<const ScalarT*[]>(static_cast<size_t>(size_));
      yp_ext_ = std::make_unique<const ScalarT*[]>(static_cast<size_t>(size_));
      f_ext_  = std::make_unique<ScalarT*[]>(static_cast<size_t>(size_));

      connection_nodes_ = std::make_unique<IdxT[]>(static_cast<size_t>(size_));

      tag_.resize(static_cast<size_t>(size_));

      if (!allocated_)
      {
        allocateVectors(size_);
      }

      allocated_ = true;
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

      f_.setDataUpdated();

      return evaluateExternalResidual();
    }

    /**
     * @brief Evaluate all residuals for the component's internal variables,
     * writing them through `f_int_`.
     *
     * @return An error code, or 0 if successful.
     */
    virtual int evaluateInternalResidual() = 0;

    /**
     * @brief Evaluate all residual contributions for the component's external variables,
     * accumulating them through `f_ext_`.
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

    /**
     * @brief Check whether the component has already been allocated.
     *
     * @return true if allocate() has previously completed, false otherwise.
     */
    bool isAllocated() const
    {
      return allocated_;
    }

  protected:
    /**
     * @brief Allocate state and residual storage owned by this component.
     *
     * Most components do not need state and residual storages. The most notable exception
     * is currently the system, so a separate flag is provided for the system.
     * Systems still can't directly access \ref y_, \ref yp_, and \ref f_, so they need
     * their corresponding \ref y_int_, \ref yp_int_, and \ref f_int_ set, since there isn't
     * another system above them to set it.
     *
     * @todo This is a weird exception specifically for systems - and in a hierarchical setting
     * will only be needed by the *topmost* system - subsystems shouldn't allocate and should have their
     * internal pointers set by the system above them. Ideally we can remove this exception by having
     * the integrator allocate these buffers instead of the system and set the internal pointers for the
     * topmost system.
     */
    void allocateVectors(IdxT n, bool system = false)
    {
      abs_tol_.resize(n);

      if (system)
      {
        y_.resize(n);
        yp_.resize(n);
        f_.resize(n);

        y_int_  = y_.getData();
        yp_int_ = yp_.getData();
        f_int_  = f_.getData();
      }
    }

    /// Number of external variables in this component - ones which are referenced but not owned by this component.
    size_t                  n_extern_;
    /// Number of internal variables in this component - ones which are only referenced by this component.
    size_t                  n_intern_;
    /**
     * @brief A set of variable indices which correspond to the external variables. Variables indices not in this set are internal.
     *
     * @invariant Must have a size of \ref n_extern_. Each element must be in the range [0, \ref size_ - 1]. Not currently verified anywhere.
     */
    std::set<IdxT>          extern_indices_;
    /**
     * @brief A map from local variable indices to system (global) variable indices. Used for Jacobian construction in
     * \ref PowerElectronicsModel::evaluateJacobian().
     * @note If a variable does not map to a corresponding variable in the system (such as with reference nodes), a special
     * sentinel value of \ref INVALID_INDEX is used. During Jacobian construction, such rows and columns will be pruned.
     */
    std::unique_ptr<IdxT[]> connection_nodes_;

  protected:
    /// The number of variables in this component. Should be equal to \ref n_extern_ plus \ref n_intern_. \see size()
    IdxT size_{0};
    /// The number of nonzero elements in this component's Jacobian. \see nnz()
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

    /**
     * An array of (input) pointers to state values for external variables.
     * \note The size of this array is equal to \ref size_, allowing you to index it with the index
     * of the variable in question (i.e. consisten with \ref extern_indices_). Therefore, accessing
     * and dereferencing the pointer in an internal variable index is undefined behavior.
     * \see setExternalConnectionNodes()
     */
    std::unique_ptr<const ScalarT*[]> y_ext_;
    /**
     * An array of (input) pointers to derivative values for external variables.
     * \note The size of this array is equal to \ref size_, allowing you to index it with the index
     * of the variable in question (i.e. consisten with \ref extern_indices_). Therefore, accessing
     * and dereferencing the pointer in an internal variable index is undefined behavior.
     * \see setExternalConnectionNodes()
     */
    std::unique_ptr<const ScalarT*[]> yp_ext_;
    /**
     * An array of (output) pointers to residuals for external variables.
     * \note The size of this array is equal to \ref size_, allowing you to index it with the index
     * of the variable in question (i.e. consisten with \ref extern_indices_). Therefore, accessing
     * and dereferencing the pointer in an internal variable index is undefined behavior.
     * \see setExternalConnectionNodes()
     */
    std::unique_ptr<ScalarT*[]>       f_ext_;

    std::vector<bool> tag_;
    VectorT           abs_tol_;

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

  private:
    /**
     * The internal buffer for state for the component. For most components, it will be empty and shouldn't be accessed.
     * Instead use \ref y_int_ for an internal variable or \ref y_ext_ for an external variable, respectively.
     * For components which want an internal buffer (such as a system), make sure that \ref y_int_ points here.
     * \see allocateVectors()
     */
    VectorT y_;
    /**
     * The internal buffer for derivatives for the component. For most components, it will be empty and shouldn't be accessed.
     * Instead use \ref yp_int_ for an internal variable or \ref yp_ext_ for an external variable, respectively.
     * For components which want an internal buffer (such as a system), make sure that \ref yp_int_ points here.
     * \see allocateVectors()
     */
    VectorT yp_;
    /**
     * The internal buffer for state for the component. For most components, it will be empty and shouldn't be accessed.
     * Instead use \ref f_int_ for an internal variable or \ref f_ext_ for an external variable, respectively.
     * For components which want an internal buffer (such as a system), make sure that \ref f_int_ points here.
     * \see allocateVectors()
     */
    VectorT f_;
  };

} // namespace GridKit
