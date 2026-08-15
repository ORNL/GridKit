

#include "InductionMotor.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant InductionMotor model
   *
   * Calls default ModelEvaluatorImpl constructor.
   * @todo create a test case utilizing the component.
   * @todo create a unit test to check correctness of component
   *
   * @tparam ScalarT - data type for scalar variables in the model
   * @tparam IdxT - integer index type for the model
   *
   * @param[in] id - unique identifier for the component
   * @param[in] Lls - stator leakage inductance
   */

  template <class ScalarT, typename IdxT>
  InductionMotor<ScalarT, IdxT>::InductionMotor(IdxT id, RealT Lls, RealT Rs, RealT Llr, RealT Rr, RealT Lms, RealT RJ, RealT P)
    : Lls_(Lls),
      Rs_(Rs),
      Llr_(Llr),
      Rr_(Rr),
      Lms_(Lms),
      RJ_(RJ),
      P_(P)
  {
    size_           = 10;
    n_intern_       = 5;
    n_extern_       = 5;
    extern_indices_ = {0, 1, 2, 3, 4};
    idc_            = id;
  }

  template <class ScalarT, typename IdxT>
  InductionMotor<ScalarT, IdxT>::~InductionMotor()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::tagDifferentiable()
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
  int InductionMotor<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
  {
    abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
    return 0;
  }

  /**
   * @brief Contributes to the resisdual
   *
   */
  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::evaluateInternalResidual()
  {
    f_int_[0] = (1.0 / 3.0) * (2.0 * *y_ext_[0] - *y_ext_[1] - *y_ext_[2]) - Rs_ * y_int_[0] - (Lls_ + Lms_) * yp_int_[0] - Lms_ * yp_int_[1];
    f_int_[1] = (1.0 / std::sqrt(3.0)) * (-*y_ext_[1] + *y_ext_[2]) - Rs_ * y_int_[1] - (Lls_ + Lms_) * yp_int_[1] - Lms_ * yp_int_[0];
    f_int_[2] = (*y_ext_[0] + *y_ext_[1] + *y_ext_[2]) / 3.0 - Rs_ * y_int_[2] - Lls_ * yp_int_[7];
    f_int_[3] = Rr_ * y_int_[3] + (Llr_ + Lms_) * yp_int_[3] + Lms_ * yp_int_[0] - (P_ / 2.0) * *y_ext_[3] * ((Llr_ + Lms_) * y_int_[4] + Lms_ * y_int_[1]);
    f_int_[4] = Rr_ * y_int_[4] + (Llr_ + Lms_) * yp_int_[4] + Lms_ * yp_int_[1] + (P_ / 2.0) * *y_ext_[3] * ((Llr_ + Lms_) * y_int_[3] + Lms_ * y_int_[0]);
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::evaluateExternalResidual()
  {
    *f_ext_[0] += y_int_[0] + y_int_[2];
    *f_ext_[1] += (-1.0 / 2.0) * y_int_[0] - (std::sqrt(3.0) / 2.0) * y_int_[1] + y_int_[2];
    *f_ext_[2] += (-1.0 / 2.0) * y_int_[0] + (std::sqrt(3.0) / 2.0) * y_int_[1] + y_int_[2];
    *f_ext_[3] += RJ_ * *yp_ext_[3] - (3.0 / 4.0) * P_ * Lms_ * (y_int_[0] * y_int_[4] - y_int_[1] * y_int_[3]);
    *f_ext_[4] += *yp_ext_[4] - *y_ext_[3];
    return 0;
  }

  /**
   * @brief Compute component Jacobian
   *
   * @todo need to implement
   *
   * @tparam ScalarT
   * @tparam IdxT
   * @return int
   */
  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::evaluateJacobian()
  {

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int InductionMotor<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  bool InductionMotor<ScalarT, IdxT>::isCloneable() const
  {
    return true;
  }

  template <class ScalarT, typename IdxT>
  CircuitComponent<ScalarT, IdxT>* InductionMotor<ScalarT, IdxT>::clone() const
  {
    return new InductionMotor<ScalarT, IdxT>(*this);
  }

  // Available template instantiations
  template class InductionMotor<double, long int>;
  template class InductionMotor<double, size_t>;
  template class InductionMotor<DependencyTracking::Variable, long int>;
  template class InductionMotor<DependencyTracking::Variable, size_t>;

} // namespace GridKit
