
#define _USE_MATH_DEFINES
#include "SynchronousMachine.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace GridKit
{

  /*!
   * @brief Constructor for a constant SynchronousMachine model
   *
   * Calls default ModelEvaluatorImpl constructor.
   * @todo This model's equations are not finished
   * @todo needs to be tested for correctness
   *
   * @tparam ScalarT - floating point type for the model
   * @tparam IdxT - integer index type for the model
   *
   * @param[in] id - unique identifier for the component
   * @param[in] Lls - stator leakage inductance
   * @param[in] Llkq - tuple of damper leakage reactances
   * @param[in] Llfd - field leakage reactance
   * @param[in] Llkd - damper leakage reactance
   * @param[in] Lmq - quadrature axis magnetizing reactance
   * @param[in] Lmd - direct axis magnetizing reactance
   * @param[in] Rs - stator resistance
   * @param[in] Rkq - tuple of damper resistances
   * @param[in] Rfd - field resistance
   * @param[in] Rkd - damper resistance
   * @param[in] RJ - rotor moment of inertia
   * @param[in] P - number of poles
   * @param[in] mub - rated frequency
   */

  template <class ScalarT, typename IdxT>
  SynchronousMachine<ScalarT, IdxT>::SynchronousMachine(IdxT id, RealT Lls, std::tuple<RealT, RealT> Llkq, RealT Llfd, RealT Llkd, RealT Lmq, RealT Lmd, RealT Rs, std::tuple<RealT, RealT> Rkq, RealT Rfd, RealT Rkd, RealT RJ, RealT P, RealT mub)
    : Lls_(Lls),
      Llkq_(Llkq),
      Llfd_(Llfd),
      Llkd_(Llkd),
      Lmq_(Lmq),
      Lmd_(Lmd),
      Rs_(Rs),
      Rkq_(Rkq),
      Rfd_(Rfd),
      Rkd_(Rkd),
      RJ_(RJ),
      P_(P),
      mub_(mub)
  {
    size_           = 13;
    n_intern_       = 6;
    n_extern_       = 7;
    extern_indices_ = {0, 1, 2, 3, 4};
    idc_            = id;
  }

  template <class ScalarT, typename IdxT>
  SynchronousMachine<ScalarT, IdxT>::~SynchronousMachine()
  {
  }

  /**
   * Initialization of the grid model
   */
  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::initialize()
  {
    return 0;
  }

  /*
   * \brief Identify differential variables
   */
  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::tagDifferentiable()
  {
    return 0;
  }

  /**
   * @brief Compute the resisdual of the component.
   *
   * @todo not finished
   */
  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::evaluateInternalResidual()
  {
    ScalarT                  rkq1  = static_cast<ScalarT>(std::get<0>(Rkq_));
    [[maybe_unused]] ScalarT rkq2  = static_cast<ScalarT>(std::get<1>(Rkq_));
    ScalarT                  llkq1 = static_cast<ScalarT>(std::get<0>(Llkq_));
    [[maybe_unused]] ScalarT llkq2 = static_cast<ScalarT>(std::get<1>(Llkq_));

    ScalarT cos1   = std::cos((P_ / 2.0) * y_[5]);
    ScalarT sin1   = std::sin((P_ / 2.0) * y_[5]);
    ScalarT cos23m = std::cos((P_ / 2.0) * y_[5] - (2.0 / 3.0) * M_PI);
    ScalarT sin23m = std::sin((P_ / 2.0) * y_[5] - (2.0 / 3.0) * M_PI);
    ScalarT cos23p = std::cos((P_ / 2.0) * y_[5] + (2.0 / 3.0) * M_PI);
    ScalarT sin23p = std::sin((P_ / 2.0) * y_[5] + (2.0 / 3.0) * M_PI);

    f_[5] = (-2.0 / 3.0) * (y_[0] * cos1 + y_[1] * cos23m + y_[2] * cos23p) + Rs_ * y_[6] + (Lls_ + Lmq_) * yp_[6] + Lmq_ * yp_[9] + Lmq_ * yp_[10] + y_[4] * (P_ / 2.0) * ((Lls_ + Lmd_) * y_[7] + Lmd_ * y_[11] + Lmd_ * y_[12]);
    f_[6] = (-2.0 / 3.0) * (y_[0] * sin1 - y_[1] * sin23m - y_[2] * sin23p) + Rs_ * y_[7] + (Lls_ + Lmd_) * yp_[7] + Lmd_ * yp_[11] + Lmd_ * yp_[12] - y_[4] * (P_ / 2.0) * ((Lls_ + Lmq_) * y_[6] + Lmq_ * y_[9] + Lmq_ * y_[10]);
    f_[7] = (-1.0 / 3.0) * (y_[0] + y_[1] + y_[2]) + Rs_ * y_[8] + Lls_ * yp_[8];
    f_[8] = rkq1 * y_[9] + (llkq1 + Lmq_) * yp_[9] + Lmq_ * yp_[6] + Lmq_ * yp_[10];
    f_[9] = rkq1 * y_[9] + (llkq1 + Lmq_) * yp_[9] + Lmq_ * yp_[6] + Lmq_ * yp_[10];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::evaluateExternalResidual()
  {
    [[maybe_unused]] ScalarT rkq2  = static_cast<ScalarT>(std::get<1>(Rkq_));
    [[maybe_unused]] ScalarT llkq2 = static_cast<ScalarT>(std::get<1>(Llkq_));

    ScalarT cos1   = std::cos((P_ / 2.0) * y_[5]);
    ScalarT sin1   = std::sin((P_ / 2.0) * y_[5]);
    ScalarT cos23m = std::cos((P_ / 2.0) * y_[5] - (2.0 / 3.0) * M_PI);
    ScalarT sin23m = std::sin((P_ / 2.0) * y_[5] - (2.0 / 3.0) * M_PI);
    ScalarT cos23p = std::cos((P_ / 2.0) * y_[5] + (2.0 / 3.0) * M_PI);
    ScalarT sin23p = std::sin((P_ / 2.0) * y_[5] + (2.0 / 3.0) * M_PI);

    f_[0] = y_[6] * cos1 + y_[7] * sin1 + y_[8];
    f_[1] = y_[6] * cos23m + y_[7] * sin23m + y_[8];
    f_[2] = y_[6] * cos23p + y_[7] * sin23p + y_[8];
    f_[3] = RJ_ * yp_[4] - (3.0 / 4.0) * P_ * (Lmd_ * y_[6] * (y_[7] + y_[11] + y_[12]) - Lmq_ * y_[7] * (y_[6] + y_[9] + y_[0]));
    f_[4] = yp_[5] - y_[4];
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::evaluateJacobian()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::evaluateIntegrand()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::initializeAdjoint()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::evaluateAdjointResidual()
  {
    return 0;
  }

  template <class ScalarT, typename IdxT>
  int SynchronousMachine<ScalarT, IdxT>::evaluateAdjointIntegrand()
  {
    return 0;
  }

  // Available template instantiations
  template class SynchronousMachine<double, long int>;
  template class SynchronousMachine<double, size_t>;
  template class SynchronousMachine<DependencyTracking::Variable, long int>;
  template class SynchronousMachine<DependencyTracking::Variable, size_t>;

} // namespace GridKit
