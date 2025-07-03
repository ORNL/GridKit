/**
 * @file GenClassicalImpl.hpp
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a Classical generator model.
 *
 *
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "GenClassical.hpp"
#include <Model/PhasorDynamics/Bus/Bus.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a classical generator model.
     *
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
      size_ = 5;
      setDerivedParams();

      // Temporary, to eliminate compiler warnings
      (void) busID_;
      (void) unit_id_;
    }

    /*!
     * @brief Constructor for a pi-model branch
     *
     */
    template <class ScalarT, typename IdxT>
    GenClassical<ScalarT, IdxT>::GenClassical(bus_type* bus,
                                              int       unit_id,
                                              real_type p0,
                                              real_type q0,
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
      size_ = 5;
      setDerivedParams();
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::allocate()
    {
      auto size = static_cast<size_t>(size_);
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      fB_.resize(size);
      yB_.resize(size);
      ypB_.resize(size);
      return 0;
    }

    /**
     * Initialization of the generator model
     *
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::initialize()
    {
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT p     = static_cast<ScalarT>(p0_);
      ScalarT q     = static_cast<ScalarT>(q0_);
      ScalarT vm2   = vr * vr + vi * vi;
      ScalarT ir    = (p * vr + q * vi) / vm2;
      ScalarT ii    = (p * vi - q * vr) / vm2;
      ScalarT Er    = Ra_ * ir - Xdp_ * ii + vr;
      ScalarT Ei    = Ra_ * ii + Xdp_ * ir + vi;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega = static_cast<ScalarT>(1.0);
      ScalarT Ep    = std::sqrt(Er * Er + Ei * Ei);
      ScalarT Te    = G_ * Ep * Ep - Ep * ((G_ * vr + -B_ * vi) * std::cos(delta) + (B_ * vr + G_ * vi) * std::sin(delta));

      y_[0]      = delta;
      y_[1]      = omega;
      y_[2]      = Te;
      y_[3]      = ir;
      y_[4]      = ii;
      pmech_set_ = Te;
      ep_set_    = Ep;

      for (size_t i = 0; i < static_cast<size_t>(size_); ++i)
        yp_[i] = 0.0;

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 2;
      }
      return 0;
    }

    /**
     * @brief Residual contribution computed locally
     *
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateResidualLocally(ScalarT* y, ScalarT* yp, ScalarT* f)
    {
      // Set variable aliases for better reliability
      const ScalarT delta = y[0];
      const ScalarT omega = y[1];
      const ScalarT telec = y[2];
      const ScalarT ir_   = y[3];
      const ScalarT ii_   = y[4];
      const ScalarT pmech = pmech_set_; /* Later optionally acquire from governor */
      const ScalarT ep    = ep_set_;    /* Later optionally acquire from exciter */

      // Set derivative aliases for better reliability
      const ScalarT delta_dot = yp[0];
      const ScalarT omega_dot = yp[1];

      // GenClassical differential equations
      f[0] = delta_dot - (omega - 1.0) * (2.0 * M_PI * 60.0);
      f[1] = omega_dot - (1.0 / (2.0 * H_)) * ((pmech - D_ * (omega - 1.0)) / omega - telec);

      // GenClassical algebraic equations
      f[2] = telec - (G_ * ep * ep - ep * ((G_ * vr_ + -B_ * vi_) * std::cos(delta) + (B_ * vr_ + G_ * vi_) * std::sin(delta)));

      f[3] = ir_ + G_ * vr_ - B_ * vi_ - ep * (G_ * std::cos(delta) - B_ * std::sin(delta));
      f[4] = ii_ + B_ * vr_ + G_ * vi_ - ep * (B_ * std::cos(delta) + G_ * std::sin(delta));

      return 0;
    }

    /**
     * \brief Residual for the generator model.
     *
     */
    template <class ScalarT, typename IdxT>
    int GenClassical<ScalarT, IdxT>::evaluateResidual()
    {
      vr_ = Vr();
      vi_ = Vi();

      evaluateResidualLocally(y_.data(), yp_.data(), f_.data());

      // GenClassical contribution to bus algebraic equations
      Ir() += ir_;
      Ii() += ii_;

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
      B_ = -Xdp_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
    }

  } // namespace PhasorDynamics
} // namespace GridKit
