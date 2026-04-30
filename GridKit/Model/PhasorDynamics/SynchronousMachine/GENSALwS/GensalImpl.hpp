#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSALwS/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSALwS/GensalData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Constructor for a GENSAL generator model with saturation
     */
    template <class ScalarT, typename IdxT>
    Gensal<ScalarT, IdxT>::Gensal(bus_type* bus, const model_data_type& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 16;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    Gensal<ScalarT, IdxT>::~Gensal()
    {
    }

    /// Helper function to extract and assign model parameters from the model's associated
    /// data structure.
    template <class ScalarT, typename IdxT>
    void Gensal<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
    {
      if (data.parameters.contains(model_data_type::Parameters::p0))
      {
        p0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::p0));
      }

      if (data.parameters.contains(model_data_type::Parameters::q0))
      {
        q0_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::q0));
      }

      if (data.parameters.contains(model_data_type::Parameters::H))
      {
        H_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::H));
      }

      if (data.parameters.contains(model_data_type::Parameters::D))
      {
        D_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::D));
      }

      if (data.parameters.contains(model_data_type::Parameters::Ra))
      {
        Ra_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Ra));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tdop))
      {
        Tdop_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tdop));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tdopp))
      {
        Tdopp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tdopp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Tqopp))
      {
        Tqopp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tqopp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xd))
      {
        Xd_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xd));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xdp))
      {
        Xdp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xdp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xdpp))
      {
        Xdpp_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xdpp));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xq))
      {
        Xq_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xq));
      }

      if (data.parameters.contains(model_data_type::Parameters::Xl))
      {
        Xl_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Xl));
      }

      if (data.parameters.contains(model_data_type::Parameters::S10))
      {
        S10_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::S10));
      }

      if (data.parameters.contains(model_data_type::Parameters::S12))
      {
        S12_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::S12));
      }

      if (data.parameters.contains(model_data_type::Parameters::mva_base))
      {
        mva_base_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::mva_base));
      }
    }

    template <class ScalarT, typename IdxT>
    const Model::VariableMonitorBase* Gensal<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

    template <class ScalarT, typename IdxT>
    void Gensal<ScalarT, IdxT>::initializeMonitor()
    {
      using Variable = typename model_data_type::MonitorableVariables;
      monitor_->set(Variable::ir, [this]
                    { return y_[12]; });
      monitor_->set(Variable::ii, [this]
                    { return y_[13]; });
      monitor_->set(Variable::p, [this]
                    { return Vr() * Ir() + Vi() * Ii(); });
      monitor_->set(Variable::q, [this]
                    { return Vi() * Ir() - Vr() * Ii(); });
      monitor_->set(Variable::delta, [this]
                    { return y_[0]; });
      monitor_->set(Variable::omega, [this]
                    { return y_[1]; });
      monitor_->set(Variable::speed, [this]
                    { return 1.0 + y_[1]; });
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int Gensal<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Gensal<ScalarT, IdxT>::allocate()
    {
      // Resize component model data
      auto size = static_cast<size_t>(size_);
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Resize bus data
      wb_.resize(2);
      h_.resize(2);

      // Resize signal variable data
      ws_.resize(2);
      ws_indices_.resize(2);
      ws_indices_[0] = INVALID_INDEX<IdxT>;
      ws_indices_[1] = INVALID_INDEX<IdxT>;

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Set output signals
      if (signals_.template isAssigned<GensalInternalVariables::OMEGA>())
      {
        signals_.template getSignalNode<GensalInternalVariables::OMEGA>()->set(&y_[1], &(this->getVariableIndex(1)));
      }

      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <class ScalarT, typename IdxT>
    int Gensal<ScalarT, IdxT>::verify() const
    {
      static constexpr auto PM  = GensalExternalVariables::PM;
      static constexpr auto EFD = GensalExternalVariables::EFD;

      int ret = 0;

      if (signals_.template isAttached<PM>())
      {
        if (!signals_.template isLinked<PM>())
        {
          Log::error() << "Gensal: pmech signal attached with no linked governor\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<EFD>())
      {
        if (!signals_.template isLinked<EFD>())
        {
          Log::error() << "Gensal: efd signal attached with no linked exciter\n";
          ret += 1;
        }
      }

      return ret;
    }

    /**
     * Initialization of the generator model
     *
     */
    template <class ScalarT, typename IdxT>
    int Gensal<ScalarT, IdxT>::initialize()
    {
      // Network frame terminal values
      ScalarT vr  = Vr();
      ScalarT vi  = Vi();
      ScalarT p   = static_cast<ScalarT>(p0_) * mva_system_base_ / mva_base_;
      ScalarT q   = static_cast<ScalarT>(q0_) * mva_system_base_ / mva_base_;
      ScalarT vm2 = vr * vr + vi * vi;
      ScalarT ir  = (p * vr + q * vi) / vm2;
      ScalarT ii  = (p * vi - q * vr) / vm2;

      ScalarT Er    = vr + Ra_ * ir - Xq_ * ii;
      ScalarT Ei    = vi + Ra_ * ii + Xq_ * ir;
      ScalarT delta = std::atan2(Ei, Er);
      ScalarT omega(0.0);

      ScalarT id      = ir * std::sin(delta) - ii * std::cos(delta);
      ScalarT iq      = ir * std::cos(delta) + ii * std::sin(delta);
      ScalarT psiqpp  = -Xq2_ * iq;
      ScalarT vd      = -psiqpp * (ONE<RealT> + omega);
      ScalarT vq      = vr * std::cos(delta) + vi * std::sin(delta) + id * Xdpp_ + iq * Ra_;
      ScalarT psidpp  = vq / (ONE<RealT> + omega);
      ScalarT psidp   = psidpp - (Xdpp_ - Xl_) * id;
      ScalarT Eqp     = psidp + Xd2_ * id;
      ScalarT Eqp_sat = Eqp - SA_;
      ScalarT ksat    = SB_ * Eqp_sat * Eqp_sat * Math::sigmoid(Eqp_sat);
      ScalarT Te      = (psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id;

      y_[0]  = delta;
      y_[1]  = omega;
      y_[2]  = Eqp;
      y_[3]  = psidp;
      y_[4]  = psiqpp;
      y_[5]  = psidpp;
      y_[6]  = ksat;
      y_[7]  = vd;
      y_[8]  = vq;
      y_[9]  = Te;
      y_[10] = id;
      y_[11] = iq;
      y_[12] = ir;
      y_[13] = ii;
      y_[14] = G_ * (vd * std::sin(delta) + vq * std::cos(delta))
               - B_ * (vd * -std::cos(delta) + vq * std::sin(delta));
      y_[15] = B_ * (vd * std::sin(delta) + vq * std::cos(delta))
               + G_ * (vd * -std::cos(delta) + vq * std::sin(delta));

      pmech_set_ = Te;
      if (signals_.template isAttached<GensalExternalVariables::PM>())
      {
        signals_.template writeExternalVariable<GensalExternalVariables::PM>(Te);
      }

      efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + ksat;
      if (signals_.template isAttached<GensalExternalVariables::EFD>())
      {
        signals_.template writeExternalVariable<GensalExternalVariables::EFD>(efd_set_);
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
    int Gensal<ScalarT, IdxT>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 5;
      }
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Gensal<ScalarT, IdxT>::evaluateInternalResidual(
        ScalarT* y,
        ScalarT* yp,
        ScalarT* wb,
        ScalarT* ws,
        ScalarT* f)
    {
      /* Read variables */
      ScalarT delta  = y[0];
      ScalarT omega  = y[1];
      ScalarT Eqp    = y[2];
      ScalarT psidp  = y[3];
      ScalarT psiqpp = y[4];
      ScalarT psidpp = y[5];
      ScalarT ksat   = y[6];
      ScalarT vd     = y[7];
      ScalarT vq     = y[8];
      ScalarT telec  = y[9];
      ScalarT id     = y[10];
      ScalarT iq     = y[11];
      ScalarT ir     = y[12];
      ScalarT ii     = y[13];
      ScalarT inr    = y[14];
      ScalarT ini    = y[15];

      /* Read derivatives */
      ScalarT delta_dot  = yp[0];
      ScalarT omega_dot  = yp[1];
      ScalarT Eqp_dot    = yp[2];
      ScalarT psidp_dot  = yp[3];
      ScalarT psiqpp_dot = yp[4];

      // Set coupling variable aliases
      ScalarT vr = wb[0];
      ScalarT vi = wb[1];

      // Set signal variable aliases
      ScalarT pmech = ws[0];
      ScalarT efd   = ws[1];

      /* 5 Gensal differential equations */
      // TODO: Replace hard-coded 60 Hz with the model's resolved frequency base
      // once frequency-base ownership/propagation is designed consistently.
      f[0] = delta_dot - omega * (TWO<RealT> * M_PI * 60.0);
      f[1] = omega_dot - (ONE<RealT> / (TWO<RealT> * H_)) * ((pmech - D_ * omega) / (ONE<RealT> + omega) - telec);
      f[2] = Eqp_dot - (ONE<RealT> / Tdop_) * (efd - (Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + ksat));
      f[3] = psidp_dot - (ONE<RealT> / Tdopp_) * (Eqp - psidp - Xd2_ * id);
      f[4] = psiqpp_dot - (ONE<RealT> / Tqopp_) * (-psiqpp - Xq2_ * iq);

      /* 9 Gensal algebraic equations */
      f[5]            = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
      ScalarT Eqp_sat = Eqp - SA_;
      f[6]            = ksat - SB_ * Eqp_sat * Eqp_sat * Math::sigmoid(Eqp_sat);
      f[7]            = vd + psiqpp * (ONE<RealT> + omega);
      f[8]            = vq - psidpp * (ONE<RealT> + omega);
      f[9]            = telec - ((psidpp - id * Xdpp_) * iq - (psiqpp - iq * Xdpp_) * id);
      f[10]           = id - (ir * std::sin(delta) - ii * std::cos(delta));
      f[11]           = iq - (ir * std::cos(delta) + ii * std::sin(delta));
      f[12]           = ir + G_ * vr - B_ * vi - inr;
      f[13]           = ii + B_ * vr + G_ * vi - ini;

      /* 2 Gensal current source definitions */
      f[14] = inr - (G_ * (std::sin(delta) * vd + std::cos(delta) * vq) - B_ * (-std::cos(delta) * vd + std::sin(delta) * vq));
      f[15] = ini - (B_ * (std::sin(delta) * vd + std::cos(delta) * vq) + G_ * (-std::cos(delta) * vd + std::sin(delta) * vq));

      return 0;
    }

    /**
     * @brief Bus residual
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Gensal<ScalarT, IdxT>::evaluateBusResidual(
        ScalarT*                  y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      ScalarT inr = y[14];
      ScalarT ini = y[15];
      ScalarT vr  = wb[0];
      ScalarT vi  = wb[1];

      h[0] = (inr - vr * G_ + vi * B_) * mva_base_ / mva_system_base_;
      h[1] = (ini - vr * B_ - vi * G_) * mva_base_ / mva_system_base_;

      return 0;
    }

    /**
     * \brief Residual evaluation and contribution to the connected bus
     *
     */
    template <class ScalarT, typename IdxT>
    int Gensal<ScalarT, IdxT>::evaluateResidual()
    {
      // Mechanical Power
      ws_[0] = pmech_set_;
      if (signals_.template isAttached<GensalExternalVariables::PM>())
      {
        ws_[0]         = signals_.template readExternalVariable<GensalExternalVariables::PM>();
        ws_indices_[0] = signals_.template readExternalVariableIndex<GensalExternalVariables::PM>();
      }

      // Exciter Efield
      ws_[1] = efd_set_;
      if (signals_.template isAttached<GensalExternalVariables::EFD>())
      {
        ws_[1]         = signals_.template readExternalVariable<GensalExternalVariables::EFD>();
        ws_indices_[1] = signals_.template readExternalVariableIndex<GensalExternalVariables::EFD>();
      }

      // Bus voltages
      wb_[0] = Vr();
      wb_[1] = Vi();

      // Residual evaluation
      evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
      evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());

      // Gensal contribution to bus algebraic equations
      Ir() += h_[0];
      Ii() += h_[1];

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
    ScalarT Gensal<ScalarT, IdxT>::getSpeed()
    {
      return y_[1];
    }

    template <class ScalarT, typename IdxT>
    ScalarT Gensal<ScalarT, IdxT>::getTorque()
    {
      return y_[9];
    }

    template <class ScalarT, typename IdxT>
    void Gensal<ScalarT, IdxT>::setDerivedParams()
    {
      SA_ = 0;
      SB_ = 0;
      if (S12_ != 0)
      {
        RealT s112 = std::sqrt(S10_ / S12_);

        SA_ = (1.2 * s112 + ONE<RealT>) / (s112 + ONE<RealT>);
        SB_ = (1.2 * s112 - ONE<RealT>) / (s112 - ONE<RealT>);
        if (SB_ < SA_)
          SA_ = SB_;
        SB_ = S12_ / ((SA_ - 1.2) * (SA_ - 1.2));
      }
      Xd1_ = Xd_ - Xdp_;
      Xd2_ = Xdp_ - Xl_;
      Xd3_ = (Xdp_ - Xdpp_) / (Xd2_ * Xd2_);
      Xd4_ = (Xdp_ - Xdpp_) / Xd2_;
      Xd5_ = (Xdpp_ - Xl_) / Xd2_;
      Xq2_ = Xq_ - Xdpp_;
      G_   = Ra_ / (Ra_ * Ra_ + Xdpp_ * Xdpp_);
      B_   = -Xdpp_ / (Ra_ * Ra_ + Xdpp_ * Xdpp_);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
