#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "Genrou.hpp"
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp> // <- TODO: Temporary, to be removed.
#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a GENROU generator model with saturation
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
        bus_id_(0),
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
      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus,
                                  IdxT      unit_id,
                                  real_type p0,
                                  real_type q0,
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
        bus_id_(0),
        unit_id_(unit_id),
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
      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, const model_data_type& data)
      : bus_(bus),
        unit_id_(1)
    {
      initializeParameters(data);

      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, signal_type* omega, signal_type* pmech, const model_data_type& data)
      : bus_(bus),
        unit_id_(1)
    {
      signals_.template attachSignalNode<GenrouExternalVariables::PM>(pmech);
      signals_.template assignSignalNode<GenrouInternalVariables::OMEGA>(omega);
      initializeParameters(data);

      size_ = 19;
      setDerivedParams();
    }

    /**
     * @brief Constructor for a GENROU generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Genrou<ScalarT, IdxT>::Genrou(bus_type* bus, signal_type* omega, signal_type* pmech, signal_type* efd, const model_data_type& data)
      : bus_(bus),
        unit_id_(1)
    {
      signals_.template attachSignalNode<GenrouExternalVariables::PM>(pmech);
      signals_.template assignSignalNode<GenrouInternalVariables::OMEGA>(omega);
      signals_.template attachSignalNode<GenrouExternalVariables::EFD>(efd);
      initializeParameters(data);

      size_ = 19;
      setDerivedParams();
    }

    /// Helper function to extract and assign model parameters from the model's associated
    /// data structure.
    template <class ScalarT, typename IdxT>
    void Genrou<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
    {
      if (data.parameters.contains(model_data_type::Parameters::p0))
      {
        p0_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::p0));
      }

      if (data.parameters.contains(model_data_type::Parameters::q0))
      {
        q0_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::q0));
      }

      if (data.parameters.contains(model_data_type::Parameters::H))
      {
        H_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::H));
      }

      if (data.parameters.contains(model_data_type::Parameters::D))
      {
        D_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::D));
      }

      if (data.parameters.contains(model_data_type::Parameters::Ra))
      {
        Ra_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Ra));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tdop))
      {
        Tdop_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Tdop));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tdopp))
      {
        Tdopp_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Tdopp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tqopp))
      {
        Tqopp_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Tqopp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tqop))
      {
        Tqop_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Tqop));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xd))
      {
        Xd_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xd));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xdp))
      {
        Xdp_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xdp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xdpp))
      {
        Xdpp_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xdpp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xq))
      {
        Xq_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xq));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xqp))
      {
        Xqp_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xqp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xqpp))
      {
        Xqpp_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xqpp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xl))
      {
        Xl_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::Xl));
      }

      if (data.parameters.contains(model_data_type::Parameters::S10))
      {
        S10_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::S10));
      }

      if (data.parameters.contains(model_data_type::Parameters::S12))
      {
        S12_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::S12));
      }

      if (data.ports.contains(model_data_type::Ports::bus))
      {
        bus_id_ = data.ports.at(model_data_type::Ports::bus);
      }
    }

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

      if (signals_.template isAssigned<GenrouInternalVariables::OMEGA>())
      {
        signals_.template getSignalNode<GenrouInternalVariables::OMEGA>()->set(&y_[1]);
      }

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
      ScalarT vr    = Vr();
      ScalarT vi    = Vi();
      ScalarT p     = static_cast<ScalarT>(p0_);
      ScalarT q     = static_cast<ScalarT>(q0_);
      ScalarT vm2   = vr * vr + vi * vi;
      ScalarT Er    = vr + (Ra_ * p * vr + Ra_ * q * vi - Xq_ * p * vi + Xq_ * q * vr) / vm2;
      ScalarT Ei    = vi + (Ra_ * p * vi - Ra_ * q * vr + Xq_ * p * vr + Xq_ * q * vi) / vm2;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega(0.0);
      ScalarT ir     = (p * vr + q * vi) / vm2;
      ScalarT ii     = (p * vi - q * vr) / vm2;
      ScalarT id     = ir * std::sin(delta) - ii * std::cos(delta);
      ScalarT iq     = ir * std::cos(delta) + ii * std::sin(delta);
      ScalarT vd     = vr * std::sin(delta) - vi * std::cos(delta) + id * Ra_ - iq * Xqpp_;
      ScalarT vq     = vr * std::cos(delta) + vi * std::sin(delta) + id * Xqpp_ - iq * Ra_;
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

      y_[0] = delta; //= 0.55399038;
      y_[1] = omega; // = 0;
      y_[2] = Eqp;   // = 0.995472581;
      y_[3] = psidp; // = 0.971299567;
      y_[4] = psiqp; // = 0.306880069;
      y_[5] = Edp;   // = 0;

      y_[6] = psiqpp = -psiqp * Xq4_ - Edp * Xq5_;
      y_[7] = psidpp = psidp * Xd4_ + Eqp * Xd5_;
      y_[8] = psipp = std::sqrt(psiqpp * psiqpp + psidpp * psidpp);
      y_[9] = ksat = SB_ * ((psipp - SA_) * (psipp - SA_));
      y_[10] = vd = -psiqpp * (1 + omega);
      y_[11] = vq = psidpp * (1 + omega);
      y_[12] = Te = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;
      y_[13]      = id;
      y_[14]      = iq;
      y_[15]      = ir;
      y_[16]      = ii;
      // y_[17] = efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat;
      y_[17]      = G_ * (vd * std::sin(delta) + vq * std::cos(delta))
               - B_ * (vd * -std::cos(delta) + vq * std::sin(delta)); /* inort, real */
      y_[18] = B_ * (vd * std::sin(delta) + vq * std::cos(delta))
               + G_ * (vd * -std::cos(delta) + vq * std::sin(delta)); /* inort, imag */

      // Set Setpoint mechanical power, which may or may not be used
      pmech_set_ = Te;
      if (signals_.template isAttached<GenrouExternalVariables::PM>())
      {
        signals_.template writeExternalVariable<GenrouExternalVariables::PM>(Te);
      }

      // Set Efield signal (may or may not exist)
      efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat;
      if (signals_.template isAttached<GenrouExternalVariables::EFD>())
      {
        signals_.template writeExternalVariable<GenrouExternalVariables::EFD>(efd_set_);
      }

      for (IdxT i = 0; i < size_; ++i)
      {
        yp_[static_cast<size_t>(i)] = 0.0;
      }

      return 0;
    }

    /**
     * \brief Identify differential variables.
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

    /**
     * @brief Residual contribution computed locally
     *
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::evaluateResidualLocally(ScalarT* y, ScalarT* yp, ScalarT* f)
    {
      /* Read variables */
      ScalarT delta  = y[0];
      ScalarT omega  = y[1];
      ScalarT Eqp    = y[2];
      ScalarT psidp  = y[3];
      ScalarT psiqp  = y[4];
      ScalarT Edp    = y[5];
      ScalarT psiqpp = y[6];
      ScalarT psidpp = y[7];
      ScalarT psipp  = y[8];
      ScalarT ksat   = y[9];
      ScalarT vd     = y[10];
      ScalarT vq     = y[11];
      ScalarT telec  = y[12];
      ScalarT id     = y[13];
      ScalarT iq     = y[14];
      ScalarT ir     = y[15];
      ScalarT ii     = y[16];
      ScalarT inr    = y[17];
      ScalarT ini    = y[18];

      /* Read derivatives */
      ScalarT delta_dot = yp[0];
      ScalarT omega_dot = yp[1];
      ScalarT Eqp_dot   = yp[2];
      ScalarT psidp_dot = yp[3];
      ScalarT psiqp_dot = yp[4];
      ScalarT Edp_dot   = yp[5];

      /* 6 Genrou differential equations */
      f[0] = delta_dot - omega * (2.0 * M_PI * 60.0);
      f[1] = omega_dot - (1.0 / (2.0 * H_)) * ((pmech_ - D_ * omega) / (1.0 + omega) - telec);
      f[2] = Eqp_dot - (1.0 / Tdop_) * (efd - (Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + psidpp * ksat));
      f[3] = psidp_dot - (1.0 / Tdopp_) * (Eqp - psidp - Xd2_ * id);
      f[4] = psiqp_dot - (1.0 / Tqopp_) * (Edp - psiqp + Xq2_ * iq);
      f[5] = Edp_dot - (1.0 / Tqop_) * (-Edp + Xqd_ * psiqpp * ksat + Xq1_ * (iq - Xq3_ * (Edp + iq * Xq2_ - psiqp)));

      /* 11 Genrou algebraic equations */
      f[6]  = psiqpp - (-psiqp * Xq4_ - Edp * Xq5_);
      f[7]  = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
      f[8]  = psipp - std::sqrt((psidpp * psidpp) + (psiqpp * psiqpp));
      f[9]  = ksat - SB_ * ((psipp - SA_) * (psipp - SA_));
      f[10] = vd + psiqpp * (1.0 + omega);
      f[11] = vq - psidpp * (1.0 + omega);
      f[12] = telec - ((psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id);
      f[13] = id - (ir * std::sin(delta) - ii * std::cos(delta));
      f[14] = iq - (ir * std::cos(delta) + ii * std::sin(delta));
      f[15] = ir + G_ * vr - B_ * vi - inr;
      f[16] = ii + B_ * vr + G_ * vi - ini;

      /* 2 Genrou control inputs are set to constant for this example */
      // f[17] = efd - efd_set_;

      /* 2 Genrou current source definitions */
      f[17] = inr - (G_ * (std::sin(delta) * vd + std::cos(delta) * vq) - B_ * (-std::cos(delta) * vd + std::sin(delta) * vq));
      f[18] = ini - (B_ * (std::sin(delta) * vd + std::cos(delta) * vq) + G_ * (-std::cos(delta) * vd + std::sin(delta) * vq));

      return 0;
    }

    /**
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int Genrou<ScalarT, IdxT>::evaluateResidual()
    {
      // Mechanmical Power
      ScalarT pmech;
      if (signals_.template isAttached<GenrouExternalVariables::PM>())
      {
        pmech = signals_.template readExternalVariable<GenrouExternalVariables::PM>();
      }
      else
      {
        pmech = pmech_set_;
      }

      // Exciter Efield
      ScalarT efd;
      if (signals_.template isAttached<GenrouExternalVariables::EFD>())
      {
        efd = signals_.template readExternalVariable<GenrouExternalVariables::EFD>();
      }
      else
      {
        efd = efd_set_;
      }
      
      // Bus voltages
      vr_ = Vr();
      vi_ = Vi();

      // Residual evaluation
      evaluateResidualLocally(y_.data(), yp_.data(), f_.data());

      // Genrou contribution to bus algebraic equations
      ScalarT inr = y_[18];
      ScalarT ini = y_[19];
      Ir() += inr - Vr() * G_ + Vi() * B_;
      Ii() += ini - Vr() * B_ - Vi() * G_;
      
      return 0;

    }

    /**
     * @brief Access generator relative speed
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    ScalarT Genrou<ScalarT, IdxT>::getSpeed()
    {
      return y_[1];
    }

    template <class ScalarT, typename IdxT>
    ScalarT Genrou<ScalarT, IdxT>::getTorque()
    {
      return y_[12];
    }

    template <class ScalarT, typename IdxT>
    void Genrou<ScalarT, IdxT>::setDerivedParams()
    {
      SA_ = 0;
      SB_ = 0;
      if (S12_ != 0)
      {
        real_type s112 = std::sqrt(S10_ / S12_);

        SA_ = (1.2 * s112 + 1.0) / (s112 + 1.0);
        SB_ = (1.2 * s112 - 1.0) / (s112 - 1.0);
        if (SB_ < SA_)
          SA_ = SB_;
        SB_ = S12_ / ((SA_ - 1.2) * (SA_ - 1.2));
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
  } // namespace PhasorDynamics
} // namespace GridKit
