

#pragma once

#include <map>
#include <set>

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
    using real_type = typename Model::Evaluator<ScalarT, IdxT>::real_type;

    CircuitComponent()  = default;
    ~CircuitComponent() = default;

    void updateTime(real_type t, real_type a)
    {
      this->time_  = t;
      this->alpha_ = a;
    }

    bool hasJacobian()
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
    IdxT getNodeConnection(IdxT local_index)
    {
      return connection_nodes_.at(local_index);
    }

  public:
    virtual IdxT size()
    {
      return size_;
    }

    virtual IdxT nnz()
    {
      return nnz_;
    }

    virtual IdxT sizeQuadrature()
    {
      return size_quad_;
    }

    virtual IdxT sizeParams()
    {
      return size_opt_;
    }

    virtual void setTolerances(real_type& rel_tol, real_type& abs_tol) const
    {
      rel_tol = rel_tol_;
      abs_tol = abs_tol_;
    }

    virtual void setMaxSteps(IdxT& msa) const
    {
      msa = max_steps_;
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

    GridKit::LinearAlgebra::COO_Matrix<real_type, IdxT>& getJacobian()
    {
      return jac_;
    }

    const GridKit::LinearAlgebra::COO_Matrix<real_type, IdxT>& getJacobian() const
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

    //@todo Fix ID naming
    IdxT getIDcomponent()
    {
      return idc_;
    }

  protected:
    size_t               n_extern_;
    size_t               n_intern_;
    std::set<IdxT>       extern_indices_;
    ///@todo may want to replace the mapping of connection_nodes to Node objects instead of IdxT. Allows for container free setup
    std::map<IdxT, IdxT> connection_nodes_;

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

    GridKit::LinearAlgebra::COO_Matrix<real_type, IdxT> jac_;

    std::vector<ScalarT> param_;
    std::vector<ScalarT> param_up_;
    std::vector<ScalarT> param_lo_;

    real_type time_;
    real_type alpha_;

    real_type rel_tol_;
    real_type abs_tol_;

    IdxT max_steps_;

    IdxT idc_;
  };

} // namespace GridKit
