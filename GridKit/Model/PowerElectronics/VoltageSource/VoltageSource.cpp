

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
  VoltageSource<ScalarT, IdxT>::VoltageSource(IdxT id, RealT V)
    : V_(V)
  {
    size_           = 3;
    n_intern_       = 1;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
  }

  template <class ScalarT, typename IdxT>
  VoltageSource<ScalarT, IdxT>::~VoltageSource()
  {
  }

  /*!
   * @brief allocate method computes sparsity pattern of the Jacobian.
   */
  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::allocate()
  {
    y_.resize(static_cast<size_t>(size_));
    yp_.resize(static_cast<size_t>(size_));
    f_.resize(static_cast<size_t>(size_));

    return 0;
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
   * \brief Specify absolute tolerance
   */
  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::setAbsoluteTolerance()
  {
    return 0;
  }

  /**
   * @brief Evaluate resisdual of component
   */
  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateResidual()
  {
    // input
    f_[0] = -y_[2];
    // ouput
    f_[1] = y_[2];
    // internal
    f_[2] = y_[1] - y_[0] - V_;
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int VoltageSource<ScalarT, IdxT>::evaluateJacobian()
  {
    // Create dF/dy
    std::vector<IdxT>  rcord{0, 1, 2, 2};
    std::vector<IdxT>  ccord{2, 2, 0, 1};
    std::vector<RealT> vals{-1.0, 1.0, -1.0, 1.0};
    jac_.setValues(rcord, ccord, vals);

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
