

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

  template <class ScalarT, typename IdxT>
  VoltageSource<ScalarT, IdxT>::VoltageSource(IdxT id, RealT V, NodeT* node1, NodeT* node2)
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

  template <class ScalarT, typename IdxT>
  VoltageSource<ScalarT, IdxT>::~VoltageSource()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Evaluate resisdual of component
   */
  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateInternalResidual()
  {
    // internal
    f_int_[0] = y_[1] - y_[0] - V_;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateExternalResidual()
  {
    // input
    f_[0] = -y_int_[0];
    // ouput
    f_[1] = y_int_[0];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();

    // Create dF/dy
    std::vector<IdxT>  rcord{0, 1, 2, 2};
    std::vector<IdxT>  ccord{2, 2, 0, 1};
    std::vector<RealT> vals{-1.0, 1.0, -1.0, 1.0};
    this->setJacValues(rcord, ccord, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::allocate()
  {
    CircuitComponent<ScalarT, IdxT>::allocate();

    this->setExternalConnectionNodes(0, node1_->getNodeConnection(0));
    this->setExternalConnectionNodes(1, node2_->getNodeConnection(0));

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class VoltageSource<double, long int>;
  template class VoltageSource<double, size_t>;
  template class VoltageSource<DependencyTracking::Variable, long int>;
  template class VoltageSource<DependencyTracking::Variable, size_t>;

} // namespace GridKit
