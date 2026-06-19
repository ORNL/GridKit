

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
   * @param[in] id - unique identifier for the component
   * @param[in] Lls - stator leakage inductance
   */

  template <typename scalar_type, typename index_type>
  InductionMotor<scalar_type, index_type>::InductionMotor(IdxT id, RealT Lls, RealT Rs, RealT Llr, RealT Rr, RealT Lms, RealT RJ, RealT P)
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

  template <typename scalar_type, typename index_type>
  InductionMotor<scalar_type, index_type>::~InductionMotor()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Contributes to the resisdual
   *
   */
  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::evaluateInternalResidual()
  {
    const auto* y = y_.getData();

    f_int_[0] = (1.0 / 3.0) * (2.0 * y[0] - y[1] - y[2]) - Rs_ * y_int_[0] - (Lls_ + Lms_) * yp_int_[0] - Lms_ * yp_int_[1];
    f_int_[1] = (1.0 / std::sqrt(3.0)) * (-y[1] + y[2]) - Rs_ * y_int_[1] - (Lls_ + Lms_) * yp_int_[1] - Lms_ * yp_int_[0];
    f_int_[2] = (y[0] + y[1] + y[2]) / 3.0 - Rs_ * y_int_[2] - Lls_ * yp_int_[7];
    f_int_[3] = Rr_ * y_int_[3] + (Llr_ + Lms_) * yp_int_[8] + Lms_ * yp_int_[0] - (P_ / 2.0) * y[3] * ((Llr_ + Lms_) * y_int_[4] + Lms_ * y_int_[1]);
    f_int_[4] = Rr_ * y_int_[4] + (Llr_ + Lms_) * yp_int_[9] + Lms_ * yp_int_[1] + (P_ / 2.0) * y[3] * ((Llr_ + Lms_) * y_int_[3] + Lms_ * y_int_[0]);
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::evaluateExternalResidual()
  {
    const auto* y  = y_.getData();
    const auto* yp = yp_.getData();
    auto*       f  = f_.getData();

    f[0] = y_int_[0] + y_int_[2];
    f[1] = (-1.0 / 2.0) * y_int_[0] - (std::sqrt(3.0) / 2.0) * y_int_[1] + y_int_[2];
    f[2] = (-1.0 / 2.0) * y_int_[0] + (std::sqrt(3.0) / 2.0) * y_int_[1] + y_int_[2];
    f[3] = RJ_ * yp[3] - (3.0 / 4.0) * P_ * Lms_ * (y[5] * y_int_[4] - y_int_[1] * y_int_[3]);
    f[4] = yp[4] - y[3];
    f_.setDataUpdated();
    return 0;
  }

  /**
   * @brief Compute component Jacobian
   *
   * @todo need to implement
   *
   * @return int
   */
  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::evaluateJacobian()
  {

    return 0;
  }

  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::evaluateIntegrand()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::initializeAdjoint()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int InductionMotor<scalar_type, index_type>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class InductionMotor<double, long int>;
  template class InductionMotor<double, size_t>;
  template class InductionMotor<DependencyTracking::Variable, long int>;
  template class InductionMotor<DependencyTracking::Variable, size_t>;

} // namespace GridKit
