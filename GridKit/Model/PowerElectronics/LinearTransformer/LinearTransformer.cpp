

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
   * @brief Computes the component resisdual
   */
  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateInternalResidual()
  {
    int_f_[0] = y_[0] - R0_ * int_[0] - L0_ * intp_[0] - M_ * intp_[1];
    int_f_[1] = y_[1] - R1_ * int_[1] - M_ * intp_[0] - L1_ * intp_[1];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int LinearTransformer<ScalarT, IdxT>::evaluateExternalResidual()
  {
    f_[0] = int_[0];
    f_[1] = int_[1];
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

  // Available template instantiations
  template class LinearTransformer<double, long int>;
  template class LinearTransformer<double, size_t>;
  template class LinearTransformer<DependencyTracking::Variable, long int>;
  template class LinearTransformer<DependencyTracking::Variable, size_t>;

} // namespace GridKit
