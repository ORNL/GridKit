/**
 * @file GenClassical.cpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @brief Definition of a Classical generator model.
 *
 *
 */

#include "GenClassical.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief Constructor for a pi-model branch
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 0
     * - Number of independent variables = 0
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::GenClassical(bus_type* bus, int unit_id)
      : bus_(bus),
        busID_(0),
        unit_id_(unit_id),
        p0_(0.0),
        q0_(0.0),
        H_(3.0),
        D_(0.0),
        Ra_(0.0),
        Xdp_(0.5)
    {
      size_ = 7;
      setDerivedParams();

      // Temporary, to eliminate compiler warnings
      (void) busID_;
      (void) unit_id_;
    }

    /*!
     * @brief Constructor for a pi-model branch
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 0
     * - Number of independent variables = 0
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::GenClassical(bus_type* bus,
                                              int       unit_id,
                                              ScalarT   p0,
                                              ScalarT   q0,
                                              real_type H,
                                              real_type D,
                                              real_type Ra,
                                              real_type Xdp)

      : bus_(bus),
        busID_(0),
        unit_id_(unit_id),
        p0_(p0),
        q0_(q0),
        H_(H),
        D_(D),
        Ra_(Ra),
        Xdp_(Xdp)
    {
      size_ = 7;
      setDerivedParams();
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::allocate()
    {
      f_.resize(size_);
      y_.resize(size_);
      yp_.resize(size_);
      tag_.resize(size_);
      fB_.resize(size_);
      yB_.resize(size_);
      ypB_.resize(size_);
      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::initialize()
    {
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT p     = p0_;
      ScalarT q     = q0_;
      ScalarT vm2   = vr * vr + vi * vi;
      ScalarT ir    = (p * vr + q * vi) / vm2;
      ScalarT ii    = (p * vi - q * vr) / vm2;
      ScalarT Er    = (G_ * (ir + G_ * vr - B_ * vi) + B_ * (ii + B_ * vr + G_ * vi)) / (G_ * G_ + B_ * B_);
      ScalarT Ei    = (-B_ * (ir + G_ * vr - B_ * vi) + G_ * (ii + B_ * vr + G_ * vi)) / (G_ * G_ + B_ * B_);
      ScalarT delta = atan2(Ei, Er);
      ScalarT omega = 0;
      ScalarT Ep    = sqrt(Er * Er + Ei * Ei);
      ScalarT Te    = G_ * Ep * Ep - Ep * ((G_ * vr - B_ * vi) * cos(delta) + (B_ * vr + G_ * vi) * sin(delta));

      y_[0] = delta;
      y_[1] = omega;
      y_[2] = Te;
      y_[3] = ir;
      y_[4] = ii;
      y_[5] = pmech_set_ = Te;
      y_[6] = ep_set_ = Ep;

      for (IdxT i = 0; i < size_; ++i)
        yp_[i] = 0.0;

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::tagDifferentiable()
    {

      return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateResidual()
    {
      /* Read variables */
      ScalarT delta = y_[0];
      ScalarT omega = y_[1];
      ScalarT telec = y_[2];
      ScalarT ir    = y_[3];
      ScalarT ii    = y_[4];
      ScalarT pmech = y_[5];
      ScalarT ep    = y_[6];

      /* Read derivatives */
      ScalarT delta_dot = yp_[0];
      ScalarT omega_dot = yp_[1];

      /* 6 GenClassical differential equations */
      f_[0] = delta_dot - (omega - 1.0) * (2.0 * M_PI * 60.0);
      f_[1] = omega_dot - (1.0 / (2 * H_)) * ((pmech - D_ * (omega - 1.0)) / omega - telec);

      /* 11 GenClassical algebraic equations */
      f_[2] = telec - (1.0 / omega) * (G_ * ep * ep - ep * ((G_ * Vr() - B_ * Vi()) * cos(delta) + (B_ * Vr() + G_ * Vi()) * sin(delta)));

      f_[3] = ir + G_ * Vr() - B_ * Vi() - ep * (G_ * cos(delta) - B_ * sin(delta));
      f_[4] = ii + B_ * Vr() + G_ * Vi() - ep * (B_ * cos(delta) + G_ * sin(delta));

      f_[5] = pmech - pmech_set_;
      f_[6] = ep - ep_set_;

      Ir() += -G_ * Vr() + B_ * Vi() + ep * (G_ * cos(delta) - B_ * sin(delta));
      Ii() += -B_ * Vr() - G_ * Vi() + ep * (B_ * cos(delta) + G_ * sin(delta));

      return 0;
    }

    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    /**
     * @brief Integrand (objective) evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateIntegrand()
    {
      // std::cout << "Evaluate Integrand for GenClassical..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint initialization not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::initializeAdjoint()
    {
      // std::cout << "Initialize adjoint for GenClassical..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint residual evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      // std::cout << "Evaluate adjoint residual for GenClassical..." << std::endl;
      return 0;
    }

    /**
     * @brief Adjoint integrand (objective) evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      // std::cout << "Evaluate adjoint Integrand for GenClassical..." << std::endl;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    void GenClassical<ScalarT, IdxT>::setDerivedParams()
    {
      G_ = Ra_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
      B_ = Xdp_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
    }

    // Available template instantiations
    template class GenClassical<double, long int>;
    template class GenClassical<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit