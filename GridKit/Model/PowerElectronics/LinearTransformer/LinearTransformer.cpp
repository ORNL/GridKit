

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
   * @param[in] id - unique identifier for the component
   * @param[in] L0 - inductance 0
   * @param[in] L1 - inductance 1
   * @param[in] R0 - resistance 0
   * @param[in] R1 - resistance 1
   * @param[in] M - mutual inductance
   */

  template <typename scalar_type, typename index_type>
  LinearTransformer<scalar_type, index_type>::LinearTransformer(IdxT id, RealT L0, RealT L1, RealT R0, RealT R1, RealT M)
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

  template <typename scalar_type, typename index_type>
  LinearTransformer<scalar_type, index_type>::~LinearTransformer()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Computes the component resisdual
   */
  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::evaluateInternalResidual()
  {
    const auto* y = y_.getData();

    f_int_[0] = y[0] - R0_ * y_int_[0] - L0_ * yp_int_[0] - M_ * yp_int_[1];
    f_int_[1] = y[1] - R1_ * y_int_[1] - M_ * yp_int_[0] - L1_ * yp_int_[1];
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::evaluateExternalResidual()
  {
    auto* f = f_.getData();

    f[0] = y_int_[0];
    f[1] = y_int_[1];
    f_.setDataUpdated();
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::evaluateJacobian()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int LinearTransformer<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class LinearTransformer<double, long int>;
  template class LinearTransformer<double, size_t>;
  template class LinearTransformer<DependencyTracking::Variable, long int>;
  template class LinearTransformer<DependencyTracking::Variable, size_t>;

} // namespace GridKit
