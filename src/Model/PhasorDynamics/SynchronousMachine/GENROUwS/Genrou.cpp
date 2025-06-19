/**
 * @file Genrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a GENROU generator model.
 *
 *
 */

#define _USE_MATH_DEFINES
#include "Genrou.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/GovernorBase.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {

    /*!
     * @brief Constructor for all custom arguments
     *
     * @param bus
     * @brief Constructor for a pi-model branch
     *
     * Arguments passed to ModelEvaluatorImpl:
     * - Number of equations = 0
     * - Number of independent variables = 0
     * - Number of quadratures = 0
     * - Number of optimization parameters = 0
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, IdxT unit_id)
      : bus_(bus),
        busID_(0),
        unit_id_(unit_id),
        p0_(0.),
        q0_(0.),
        H_(3.),
        D_(0.),
        Ra_(0.),
        Tdop_(7.),
        Tdopp_(.04),
        Tqopp_(.05),
        Tqop_(.75),
        Xd_(2.1),
        Xdp_(0.2),
        Xdpp_(0.18),
        Xq_(.5),
        Xqp_(.5),
        Xqpp_(.18),
        Xl_(.15),
        S10_(0.),
        S12_(0.)
    {
      size_ = 20; // 21; 20 without Pmech
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
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus,
                                  IdxT      unit_id,
                                  gov_type* gov,
                                  ScalarT   p0,
                                  ScalarT   q0,
                                  real_type H,
                                  real_type D,
                                  real_type Ra,
                                  real_type Tdop,
                                  real_type Tdopp,
                                  real_type Tqopp,
                                  real_type Tqop,
                                  real_type Xd,
                                  real_type Xdp,
                                  real_type Xdpp,
                                  real_type Xq,
                                  real_type Xqp,
                                  real_type Xqpp,
                                  real_type Xl,
                                  real_type S10,
                                  real_type S12)
      : bus_(bus),
        busID_(0),
        unit_id_(unit_id),
        gov_(gov),
        p0_(p0),
        q0_(q0),
        H_(H),
        D_(D),
        Ra_(Ra),
        Tdop_(Tdop),
        Tdopp_(Tdopp),
        Tqopp_(Tqopp),
        Tqop_(Tqop),
        Xd_(Xd),
        Xdp_(Xdp),
        Xdpp_(Xdpp),
        Xq_(Xq),
        Xqp_(Xqp),
        Xqpp_(Xqpp),
        Xl_(Xl),
        S10_(S10),
        S12_(S12)
    {
      size_ = 20;
      setDerivedParams();
    }

    /*!
     * @brief Constructor for the GENROU generator with saturation.
     *
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, const model_data_type& data)
      : bus_(bus),
        busID_(0),
        unit_id_(data.unit_id),
        p0_(data.p0),
        q0_(data.q0),
        H_(data.H),
        D_(data.D),
        Ra_(data.Ra),
        Tdop_(data.Tdop),
        Tdopp_(data.Tdopp),
        Tqopp_(data.Tqopp),
        Tqop_(data.Tqop),
        Xd_(data.Xd),
        Xdp_(data.Xdp),
        Xdpp_(data.Xdpp),
        Xq_(data.Xq),
        Xqp_(data.Xqp),
        Xqpp_(data.Xqpp),
        Xl_(data.Xl),
        S10_(data.S10),
        S12_(data.S12)
    {
      size_ = 21;
      setDerivedParams();
    }

    // /**
    //  * @brief Destroy the Genrou
    //  *
    //  * @tparam ScalarT
    //  * @tparam IdxT
    //  */
    // template <class ScalarT, typename IdxT>
    // Genrou<ScalarT, IdxT>::~Genrou()
    // {
    //   // std::cout << "Destroy Genrou..." << std::endl;
    // }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::allocate()
    {
      f_.resize(static_cast<size_t>(size_));
      y_.resize(static_cast<size_t>(size_));
      yp_.resize(static_cast<size_t>(size_));
      tag_.resize(static_cast<size_t>(size_));
      fB_.resize(static_cast<size_t>(size_));
      yB_.resize(static_cast<size_t>(size_));
      ypB_.resize(static_cast<size_t>(size_));
      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::initialize()
    {
      /* Initialization tricks -- assuming NO saturation */
      ScalarT vr     = Vr();
      ScalarT vi     = Vi();
      ScalarT p      = p0_;
      ScalarT q      = q0_;
      ScalarT vm2    = vr * vr + vi * vi;
      ScalarT Er     = vr + (Ra_ * p * vr + Ra_ * q * vi - Xq_ * p * vi + Xq_ * q * vr) / vm2;
      ScalarT Ei     = vi + (Ra_ * p * vi - Ra_ * q * vr + Xq_ * p * vr + Xq_ * q * vi) / vm2;
      ScalarT delta  = atan2(Ei, Er);
      ScalarT omega  = 0;
      ScalarT ir     = (p * vr + q * vi) / vm2;
      ScalarT ii     = (p * vi - q * vr) / vm2;
      ScalarT id     = ir * sin(delta) - ii * cos(delta);
      ScalarT iq     = ir * cos(delta) + ii * sin(delta);
      ScalarT vd     = vr * sin(delta) - vi * cos(delta) + id * Ra_ - iq * Xqpp_;
      ScalarT vq     = vr * cos(delta) + vi * sin(delta) + id * Xqpp_ - iq * Ra_;
      ScalarT psiqpp = -vd / (1 + omega);
      ScalarT psidpp = vq / (1 + omega);
      ScalarT Te     = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      ScalarT psiqp  = -(-(Xqp_ - Xl_) * iq + psiqpp * (Xqp_ - Xl_) / (Xqpp_ - Xl_))
                      / (1 + (Xqp_ - Xqpp_) / (Xqpp_ - Xl_));
      ScalarT Edp   = psiqp - (Xqp_ - Xl_) * iq;
      ScalarT psidp = -((Xdp_ - Xl_) * id - psidpp * (Xdp_ - Xl_) / (Xdpp_ - Xl_))
                      / (1 + (Xdp_ - Xdpp_) / (Xdpp_ - Xl_));
      ScalarT Eqp = psidp + (Xdp_ - Xl_) * id;

      /* Now we have the state variables, solve for alg. variables */
      ScalarT ksat;
      ScalarT psipp;

      y_[0] = delta; // = 0.55399038;
      y_[1] = omega; // = 0;
      y_[2] = Eqp;   // = 0.995472581;
      y_[3] = psidp; // = 0.971299567;
      y_[4] = psiqp; // = 0.306880069;
      y_[5] = Edp;   // = 0;

      y_[6] = psiqpp = -psiqp * Xq4_ - Edp * Xq5_;
      y_[7] = psidpp = psidp * Xd4_ + Eqp * Xd5_;
      y_[8] = psipp = sqrt(psiqpp * psiqpp + psidpp * psidpp);
      y_[9] = ksat = SB_ * pow(psipp - SA_, 2);
      y_[10] = vd = -psiqpp * (1 + omega);
      y_[11] = vq = psidpp * (1 + omega);
      y_[12] = Te = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      y_[13]      = id;
      y_[14]      = iq;
      y_[15]      = ir;
      y_[16]      = ii;
      // y_[17] = pmech_set_ = Te;
      y_[17] = efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat;
      y_[18]            = G_ * (vd * sin(delta) + vq * cos(delta))
               - B_ * (vd * -cos(delta) + vq * sin(delta)); /* inort, real */
      y_[19] = B_ * (vd * sin(delta) + vq * cos(delta))
               + G_ * (vd * -cos(delta) + vq * sin(delta)); /* inort, imag */

      // Set Setpoint mechanical power, which may or may not be used
      pmech_set_ = Te;

      for (IdxT i = 0; i < size_; ++i)
        yp_[static_cast<size_t>(i)] = 0.0;

      return 0;
    }

    /*!
     * @brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 6;
      }
      return 0;
    }

    /*!
     * @brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::evaluateResidual()
    {
      /* Read variables */
      ScalarT delta  = y_[0];
      ScalarT omega  = y_[1];
      ScalarT Eqp    = y_[2];
      ScalarT psidp  = y_[3];
      ScalarT psiqp  = y_[4];
      ScalarT Edp    = y_[5];
      ScalarT psiqpp = y_[6];
      ScalarT psidpp = y_[7];
      ScalarT psipp  = y_[8];
      ScalarT ksat   = y_[9];
      ScalarT vd     = y_[10];
      ScalarT vq     = y_[11];
      ScalarT telec  = y_[12];
      ScalarT id     = y_[13];
      ScalarT iq     = y_[14];
      ScalarT ir     = y_[15];
      ScalarT ii     = y_[16];
      ScalarT efd    = y_[17];
      ScalarT inr    = y_[18];
      ScalarT ini    = y_[19];
      ScalarT vr     = Vr();
      ScalarT vi     = Vi();
      ScalarT pmech;
      if (gov_)
      {
        pmech = gov_->Pmech(); // ISSUE IS HERE?
      }
      else
      {
        pmech = pmech_set_;
      }

      /* Read derivatives */
      ScalarT delta_dot = yp_[0];
      ScalarT omega_dot = yp_[1];
      ScalarT Eqp_dot   = yp_[2];
      ScalarT psidp_dot = yp_[3];
      ScalarT psiqp_dot = yp_[4];
      ScalarT Edp_dot   = yp_[5];

      /* 6 Genrou differential equations */
      f_[0] = delta_dot - omega * (2 * M_PI * 60);
      f_[1] = omega_dot - (1 / (2 * H_)) * ((pmech - D_ * omega) / (1 + omega) - telec);
      f_[2] = Eqp_dot - (1 / Tdop_) * (efd - (Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat));
      f_[3] = psidp_dot - (1 / Tdopp_) * (Eqp - psidp - Xd2_ * id);
      f_[4] = psiqp_dot - (1 / Tqopp_) * (Edp - psiqp + Xq2_ * iq);
      f_[5] = Edp_dot - (1 / Tqop_) * (-Edp + Xqd_ * psiqpp * ksat + Xq1_ * (iq - Xq3_ * (Edp + iq * Xq2_ - psiqp)));

      /* 11 Genrou algebraic equations */
      f_[6]  = psiqpp - (-psiqp * Xq4_ - Edp * Xq5_);
      f_[7]  = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
      f_[8]  = psipp - sqrt(pow(psidpp, 2.0) + pow(psiqpp, 2.0));
      f_[9]  = ksat - SB_ * pow(psipp - SA_, 2.0);
      f_[10] = vd + psiqpp * (1 + omega);
      f_[11] = vq - psidpp * (1 + omega);
      f_[12] = telec - ((psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id);
      f_[13] = id - (ir * sin(delta) - ii * cos(delta));
      f_[14] = iq - (ir * cos(delta) + ii * sin(delta));
      f_[15] = ir + G_ * vr - B_ * vi - inr;
      f_[16] = ii + B_ * vr + G_ * vi - ini;

      /* 2 Genrou control inputs are set to constant for this example */
      // f_[17] = pmech - pmech_set_;
      f_[17] = efd - efd_set_;

      /* 2 Genrou current source definitions */
      f_[18] = inr - (G_ * (sin(delta) * vd + cos(delta) * vq) - B_ * (-cos(delta) * vd + sin(delta) * vq));
      f_[19] = ini - (B_ * (sin(delta) * vd + cos(delta) * vq) + G_ * (-cos(delta) * vd + sin(delta) * vq));

      /* Current balance */
      Ir() += inr - Vr() * G_ + Vi() * B_;
      Ii() += ini - Vr() * B_ - Vi() * G_;

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
    int Genrou<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Genrou..." << std::endl;
      std::cout << "Jacobian evaluation not implemented!" << std::endl;
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
    int Genrou<ScalarT, IdxT>::evaluateIntegrand()
    {
      // std::cout << "Evaluate Integrand for Genrou..." << std::endl;
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
    int Genrou<ScalarT, IdxT>::initializeAdjoint()
    {
      // std::cout << "Initialize adjoint for Genrou..." << std::endl;
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
    int Genrou<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      // std::cout << "Evaluate adjoint residual for Genrou..." << std::endl;
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
    int Genrou<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      // std::cout << "Evaluate adjoint Integrand for Genrou..." << std::endl;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    void Genrou<ScalarT, IdxT>::setDerivedParams()
    {
      SA_ = 0;
      SB_ = 0;
      if (S12_ != 0)
      {
        real_type s112 = sqrt(S10_ / S12_);

        SA_ = (1.2 * s112 + 1) / (s112 + 1);
        SB_ = (1.2 * s112 - 1) / (s112 - 1);
        if (SB_ < SA_)
          SA_ = SB_;
        SB_ = S12_ / pow(SA_ - 1.2, 2);
      }
      Xd1_ = Xd_ - Xdp_;
      Xd2_ = Xdp_ - Xl_;
      Xd3_ = (Xdp_ - Xdpp_) / (Xd2_ * Xd2_);
      Xd4_ = (Xdp_ - Xdpp_) / Xd2_;
      Xd5_ = (Xdpp_ - Xl_) / Xd2_;
      Xq1_ = Xq_ - Xqp_;
      Xq2_ = Xqp_ - Xl_;
      Xq3_ = (Xqp_ - Xqpp_) / (Xq2_ * Xq2_);
      Xq4_ = (Xqp_ - Xqpp_) / Xq2_;
      Xq5_ = (Xqpp_ - Xl_) / Xq2_;
      Xqd_ = (Xq_ - Xl_) / (Xd_ - Xl_);
      G_   = Ra_ / (Ra_ * Ra_ + Xqpp_ * Xqpp_);
      B_   = -Xqpp_ / (Ra_ * Ra_ + Xqpp_ * Xqpp_);
    }

    template <class ScalarT, typename IdxT>
    ScalarT Genrou<ScalarT, IdxT>::speed()
    {
      return y_[1];
    }

    // Available template instantiations
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
