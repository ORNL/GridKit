

#include "Capacitor.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for Capacitor
   *
   * @todo this needs to be tested on some circuit
   *
   * Calls default ModelEvaluatorImpl constructor.
   */

  template <typename scalar_type, typename index_type>
  Capacitor<scalar_type, index_type>::Capacitor(IdxT id, RealT C)
    : C_(C)
  {
    size_           = 3;
    n_intern_       = 1;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
    nnz_            = 5;
  }

  template <typename scalar_type, typename index_type>
  Capacitor<scalar_type, index_type>::~Capacitor()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Evaluate the resisdual of the Capcitor
   *
   */
  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::evaluateInternalResidual()
  {
    const auto* y = y_.getData();

    f_int_[0] = -C_ * yp_int_[0] + y[0] - y[1] - y_int_[0];
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::evaluateExternalResidual()
  {
    auto* f = f_.getData();

    // input
    f[0] = C_ * yp_int_[0];
    // output
    f[1] = -C_ * yp_int_[0];
    f_.setDataUpdated();
    return 0;
  }

  /**
   * @brief Compute the Jacobian dF/dy - a dF/dy'
   *
   * @return int
   */
  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::evaluateJacobian()
  {
    this->zeroJacMatrix();
    // Create dF/dy
    std::vector<IdxT>  rcord{0, 1, 2, 2, 2};
    std::vector<IdxT>  ccord{2, 2, 0, 1, 2};
    std::vector<RealT> vals{C_ * alpha_, -C_ * alpha_, 1.0, -1.0, -1.0 - C_ * alpha_};
    this->setJacValues(rcord, ccord, vals);

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int Capacitor<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class Capacitor<double, long int>;
  template class Capacitor<double, size_t>;
  template class Capacitor<DependencyTracking::Variable, long int>;
  template class Capacitor<DependencyTracking::Variable, size_t>;

} // namespace GridKit
