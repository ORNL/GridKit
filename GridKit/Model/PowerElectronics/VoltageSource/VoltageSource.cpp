

#include "VoltageSource.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant VoltageSource model
   *
   * Calls default ModelEvaluatorImpl constructor.
   */

  template <typename scalar_type, typename index_type>
  VoltageSource<scalar_type, index_type>::VoltageSource(IdxT id, RealT V, NodeT* node1, NodeT* node2)
    : V_(V), node1_(node1), node2_(node2)
  {
    assert(node1_->size() == 1);
    assert(node2_->size() == 1);
    size_           = 3;
    n_intern_       = 1;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
    nnz_            = 4;
  }

  template <typename scalar_type, typename index_type>
  VoltageSource<scalar_type, index_type>::~VoltageSource()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Compute the absolute tolerance for each variable in the model
   *
   * @param rel_tol The relative tolerance which can be used to pick the
   *        absolute tolerance.
   * @return int 0 if successful, non-zero otherwise.
   *
   * This represents a "noise" level close to zero for which pure relative
   * error cannot be used.
   */
  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Evaluate resisdual of component
   */
  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::evaluateInternalResidual()
  {
    // internal
    const auto* y = y_.getData();

    f_int_[0] = y[1] - y[0] - V_;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::evaluateExternalResidual()
  {
    auto* f = f_.getData();

    // input
    f[0] = -y_int_[0];
    // ouput
    f[1] = y_int_[0];
    f_.setDataUpdated();
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rcord{0, 1, 2, 2};
    std::vector<IdxT>  ccord{2, 2, 0, 1};
    std::vector<RealT> vals{-1.0, 1.0, -1.0, 1.0};
    this->setJacValues(rcord, ccord, vals);

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node1_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node2_->getNodeConnection(0));

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int VoltageSource<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class VoltageSource<double, long int>;
  template class VoltageSource<double, size_t>;
  template class VoltageSource<DependencyTracking::Variable, long int>;
  template class VoltageSource<DependencyTracking::Variable, size_t>;

} // namespace GridKit
