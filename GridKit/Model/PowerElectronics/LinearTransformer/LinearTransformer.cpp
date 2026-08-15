

#include "LinearTransformer.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a LinearTransformer model
   *
   * Calls default ModelEvaluatorImpl constructor.
   * @todo Not tested in any model yet. Should be
   * @todo Has not been tested for correctness
   *
   * @tparam ScalarT - floating point type for the model
   * @tparam IdxT - integer index type for the model
   *
   * @param[in] id - unique identifier for the component
   * @param[in] L0 - inductance 0
   * @param[in] L1 - inductance 1
   * @param[in] R0 - resistance 0
   * @param[in] R1 - resistance 1
   * @param[in] M - mutual inductance
   */

  template <class ScalarT, typename IdxT>
  LinearTransformer<ScalarT, IdxT>::LinearTransformer(IdxT id, RealT L0, RealT L1, RealT R0, RealT R1, RealT M)
    : L0_(L0),
      L1_(L1),
      R0_(R0),
      R1_(R1),
      M_(M)
  {
    size_           = 4;
    n_intern_       = 2;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
  }

  template <class ScalarT, typename IdxT>
  LinearTransformer<ScalarT, IdxT>::~LinearTransformer()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::tagDifferentiable()
  {
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
  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Computes the component resisdual
   */
  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateInternalResidual()
  {
    f_int_[0] = *y_ext_[0] - R0_ * y_int_[0] - L0_ * yp_int_[0] - M_ * yp_int_[1];
    f_int_[1] = *y_ext_[1] - R1_ * y_int_[1] - M_ * yp_int_[0] - L1_ * yp_int_[1];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateExternalResidual()
  {
    *f_ext_[0] += y_int_[0];
    *f_ext_[1] += y_int_[1];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  bool LinearTransformer<ScalarT, IdxT>::isCloneable() const
  {
    return true;
  }

  template <class ScalarT, typename IdxT>
  CircuitComponent<ScalarT, IdxT>* LinearTransformer<ScalarT, IdxT>::clone() const
  {
    return new LinearTransformer<ScalarT, IdxT>(*this);
  }

  // Available template instantiations
  template class LinearTransformer<double, long int>;
  template class LinearTransformer<double, size_t>;
  template class LinearTransformer<DependencyTracking::Variable, long int>;
  template class LinearTransformer<DependencyTracking::Variable, size_t>;

} // namespace GridKit
