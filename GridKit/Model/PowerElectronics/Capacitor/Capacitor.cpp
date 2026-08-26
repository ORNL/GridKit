

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

  template <class ScalarT, typename IdxT>
  Capacitor<ScalarT, IdxT>::Capacitor(IdxT id, RealT C)
    : C_(C)
  {
    size_           = 3;
    n_intern_       = 1;
    n_extern_       = 2;
    extern_indices_ = {0, 1};
    idc_            = id;
    nnz_            = 5;
  }

  template <class ScalarT, typename IdxT>
  Capacitor<ScalarT, IdxT>::~Capacitor()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::tagDifferentiable()
  {
    // All variables are differentials
    std::fill(tag_.begin(), tag_.end(), true);
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
  int Capacitor<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Evaluate the resisdual of the Capcitor
   *
   */
  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::evaluateInternalResidual()
  {
    f_int_[0] = -C_ * yp_int_[0] + *y_ext_[0] - *y_ext_[1] - y_int_[0];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::evaluateExternalResidual()
  {
    // input
    *f_ext_[0] += C_ * yp_int_[0];
    // output
    *f_ext_[1] += -C_ * yp_int_[0];
    return 0;
  }

  /**
   * @brief Compute the Jacobian dF/dy - a dF/dy'
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::evaluateJacobian()
  {
    this->zeroJacMatrix();
    // Create dF/dy
    std::vector<IdxT>  rcord{0, 1, 2, 2, 2};
    std::vector<IdxT>  ccord{2, 2, 0, 1, 2};
    std::vector<RealT> vals{C_ * alpha_, -C_ * alpha_, 1.0, -1.0, -1.0 - C_ * alpha_};
    this->setJacValues(rcord, ccord, vals);

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Capacitor<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  bool Capacitor<ScalarT, IdxT>::isCloneable() const
  {
    return true;
  }

  template <class ScalarT, typename IdxT>
  CircuitComponent<ScalarT, IdxT>* Capacitor<ScalarT, IdxT>::clone() const
  {
    return new Capacitor<ScalarT, IdxT>(*this);
  }

  // Available template instantiations
  template class Capacitor<double, long int>;
  template class Capacitor<double, size_t>;
  template class Capacitor<DependencyTracking::Variable, long int>;
  template class Capacitor<DependencyTracking::Variable, size_t>;

} // namespace GridKit
