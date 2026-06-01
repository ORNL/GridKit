

#pragma once

#include <set>

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
    using MatrixT    = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;
    using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;

    CircuitComponent() = default;

    ~CircuitComponent()
    {
      if (connection_nodes_ != nullptr)
      {
        delete[] connection_nodes_;
      }
    };

    CircuitComponent(const CircuitComponent& other)
      : n_extern_(other.n_extern_),
        n_intern_(other.n_intern_),
        extern_indices_(other.extern_indices_),
        size_(other.size_),
        nnz_(other.nnz_),
        size_quad_(other.size_quad_),
        size_opt_(other.size_opt_),
        y_(other.y_),
        yp_(other.yp_),
        tag_(other.tag_),
        f_(other.f_),
        g_(other.g_),
        yB_(other.yB_),
        ypB_(other.ypB_),
        fB_(other.fB_),
        gB_(other.gB_),
        jac_(other.jac_),
        param_(other.param_),
        param_up_(other.param_up_),
        param_lo_(other.param_lo_),
        time_(other.time_),
        alpha_(other.alpha_),
        rel_tol_(other.rel_tol_),
        abs_tol_(other.abs_tol_),
        max_steps_(other.max_steps_),
        idc_(other.idc_)
    {
      if (other.connection_nodes_ != nullptr)
      {
        connection_nodes_ = new IdxT[size_];

        std::copy(
            other.connection_nodes_,
            other.connection_nodes_ + size_,
            connection_nodes_);
      }
      else
      {
        connection_nodes_ = nullptr;
      }
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
     * @brief Create the mappings from local to global indices
     *
     * @param local_index
     * @param global_index
     * @return int
     */
    int setExternalConnectionNodes(IdxT local_index, IdxT global_index)
    {
      if (connection_nodes_ == nullptr)
      {
        connection_nodes_ = new IdxT[size_];
      }

      connection_nodes_[local_index] = global_index;
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
      return jac_;
    }

    const MatrixT& getJacobian() const final
    {
      return jac_;
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

    //@todo Fix ID naming
    IdxT getIDcomponent() const
    {
      return idc_;
    }

  protected:
    size_t         n_extern_;
    size_t         n_intern_;
    std::set<IdxT> extern_indices_;
    ///@todo may want to replace the mapping of connection_nodes to Node objects instead of IdxT. Allows for container free setup
    IdxT*          connection_nodes_ = nullptr;

  protected:
    IdxT size_{0};
    IdxT nnz_{0};
    IdxT size_quad_{0};
    IdxT size_opt_{0};

    std::vector<ScalarT> y_;
    std::vector<ScalarT> yp_;
    std::vector<bool>    tag_;
    std::vector<ScalarT> f_;
    std::vector<ScalarT> g_;

    std::vector<ScalarT> yB_;
    std::vector<ScalarT> ypB_;
    std::vector<ScalarT> fB_;
    std::vector<ScalarT> gB_;

    MatrixT jac_;

    std::vector<ScalarT> param_;
    std::vector<ScalarT> param_up_;
    std::vector<ScalarT> param_lo_;

    RealT time_;
    RealT alpha_;

    RealT rel_tol_;
    RealT abs_tol_;

    IdxT max_steps_;

    IdxT idc_;
  };

} // namespace GridKit
