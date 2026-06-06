#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/GensalData.hpp>
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
    template <typename scalar_type, typename index_type>
    Gensal<scalar_type, index_type>::Gensal(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 16;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Gensal<scalar_type, index_type>::~Gensal()
    {
    }

    /// Helper function to extract and assign model parameters from the model's associated
    /// data structure.
    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::p0))
      {
        p0_ = std::get<RealT>(data.parameters.at(Parameter::p0));
      }

      if (data.parameters.contains(Parameter::q0))
      {
        q0_ = std::get<RealT>(data.parameters.at(Parameter::q0));
      }

      if (data.parameters.contains(Parameter::H))
      {
        H_ = std::get<RealT>(data.parameters.at(Parameter::H));
      }

      if (data.parameters.contains(Parameter::D))
      {
        D_ = std::get<RealT>(data.parameters.at(Parameter::D));
      }

      if (data.parameters.contains(Parameter::Ra))
      {
        Ra_ = std::get<RealT>(data.parameters.at(Parameter::Ra));
      }

      if (data.parameters.contains(Parameter::Tdop))
      {
        Tdop_ = std::get<RealT>(data.parameters.at(Parameter::Tdop));
      }

      if (data.parameters.contains(Parameter::Tdopp))
      {
        Tdopp_ = std::get<RealT>(data.parameters.at(Parameter::Tdopp));
      }

      if (data.parameters.contains(Parameter::Tqopp))
      {
        Tqopp_ = std::get<RealT>(data.parameters.at(Parameter::Tqopp));
      }

      if (data.parameters.contains(Parameter::Xd))
      {
        Xd_ = std::get<RealT>(data.parameters.at(Parameter::Xd));
      }

      if (data.parameters.contains(Parameter::Xdp))
      {
        Xdp_ = std::get<RealT>(data.parameters.at(Parameter::Xdp));
      }

      if (data.parameters.contains(Parameter::Xdpp))
      {
        Xdpp_ = std::get<RealT>(data.parameters.at(Parameter::Xdpp));
      }

      if (data.parameters.contains(Parameter::Xq))
      {
        Xq_ = std::get<RealT>(data.parameters.at(Parameter::Xq));
      }

      if (data.parameters.contains(Parameter::Xl))
      {
        Xl_ = std::get<RealT>(data.parameters.at(Parameter::Xl));
      }

      if (data.parameters.contains(Parameter::S10))
      {
        S10_ = std::get<RealT>(data.parameters.at(Parameter::S10));
      }

      if (data.parameters.contains(Parameter::S12))
      {
        S12_ = std::get<RealT>(data.parameters.at(Parameter::S12));
      }

      if (data.parameters.contains(Parameter::mva))
      {
        mva_base_ = std::get<RealT>(data.parameters.at(Parameter::mva));
      }
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Gensal<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    // System base -> machine base when reading system values.
    template <typename scalar_type, typename index_type>
    scalar_type Gensal<scalar_type, index_type>::toMachineBase(ScalarT value) const
    {
      return value * va_system_base_ / va_machine_base_;
    }

    // Machine base -> system base for network and signal output.
    template <typename scalar_type, typename index_type>
    scalar_type Gensal<scalar_type, index_type>::toSystemBase(ScalarT value) const
    {
      return value / toMachineBase(static_cast<ScalarT>(ONE<RealT>));
    }

    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      // Convert monitored terminal values to system base.
      monitor_->set(Variable::ir, [this]
                    { return toSystemBase(y_[12]); });
      monitor_->set(Variable::ii, [this]
                    { return toSystemBase(y_[13]); });
      monitor_->set(Variable::p, [this]
                    { return toSystemBase(Vr() * y_[12] + Vi() * y_[13]); });
      monitor_->set(Variable::q, [this]
                    { return toSystemBase(Vi() * y_[12] - Vr() * y_[13]); });
      monitor_->set(Variable::delta, [this]
                    { return y_[0]; });
      monitor_->set(Variable::omega, [this]
                    { return y_[1]; });
      monitor_->set(Variable::speed, [this]
                    { return 1.0 + y_[1]; });
      monitor_->set(Variable::Eqp, [this]
                    { return y_[2]; });
      monitor_->set(Variable::psidp, [this]
                    { return y_[3]; });
      monitor_->set(Variable::psiqpp, [this]
                    { return y_[4]; });
      monitor_->set(Variable::psidpp, [this]
                    { return y_[5]; });
      monitor_->set(Variable::vd, [this]
                    { return y_[7]; });
      monitor_->set(Variable::vq, [this]
                    { return y_[8]; });
      monitor_->set(Variable::te, [this]
                    { return y_[9]; });
      monitor_->set(Variable::id, [this]
                    { return y_[10]; });
      monitor_->set(Variable::iq, [this]
                    { return y_[11]; });
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::allocate()
    {
      // Resize component model data
      auto size = static_cast<size_t>(size_);
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);
      abs_tol_.resize(size);
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
      if (signals_.template isAssigned<GensalOutputPorts::speed>())
      {
        signals_.template getSignalNode<GensalOutputPorts::speed>()->set(&y_[1], &(this->getVariableIndex(1)));
      }

      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::verify() const
    {
      static constexpr auto PM  = GensalInputPorts::pmech;
      static constexpr auto EFD = GensalInputPorts::efd;

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
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::initialize()
    {
      // Network frame terminal values
      ScalarT vr  = Vr();
      ScalarT vi  = Vi();
      ScalarT p   = toMachineBase(static_cast<ScalarT>(p0_));
      ScalarT q   = toMachineBase(static_cast<ScalarT>(q0_));
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

      // Convert Te to system base for governor PM signal.
      pmech_set_ = toSystemBase(Te);
      if (signals_.template isAttached<GensalInputPorts::pmech>())
      {
        signals_.template writeExternalVariable<GensalInputPorts::pmech>(pmech_set_);
      }

      efd_set_ = Eqp + Xd1_ * (id + Xd3_ * (Eqp - psidp - Xd2_ * id)) + ksat;
      if (signals_.template isAttached<GensalInputPorts::efd>())
      {
        signals_.template writeExternalVariable<GensalInputPorts::efd>(efd_set_);
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
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::tagDifferentiable()
    {
      for (IdxT i = 0; i < size_; ++i)
      {
        tag_[static_cast<size_t>(i)] = i < 5;
      }
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Gensal<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y,
        const ScalarT* yp,
        const ScalarT* wb,
        const ScalarT* ws,
        ScalarT*       f)
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
      ScalarT pmech = toMachineBase(ws[0]);
      ScalarT efd   = ws[1];

      /* 5 Gensal differential equations */
      f[0] = delta_dot - omega * (TWO<RealT> * M_PI * freq_system_base_);
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
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Gensal<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* wb,
        ScalarT*                        h)
    {
      ScalarT ir = y[12];
      ScalarT ii = y[13];

      // Convert current injection to system base for the network.
      h[0] = toSystemBase(ir);
      h[1] = toSystemBase(ii);

      return 0;
    }

    /**
     * \brief Residual evaluation and contribution to the connected bus
     *
     */
    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::evaluateResidual()
    {
      // Mechanical Power
      ws_[0] = pmech_set_;
      if (signals_.template isAttached<GensalInputPorts::pmech>())
      {
        ws_[0]         = signals_.template readExternalVariable<GensalInputPorts::pmech>();
        ws_indices_[0] = signals_.template readExternalVariableIndex<GensalInputPorts::pmech>();
      }

      // Exciter Efield
      ws_[1] = efd_set_;
      if (signals_.template isAttached<GensalInputPorts::efd>())
      {
        ws_[1]         = signals_.template readExternalVariable<GensalInputPorts::efd>();
        ws_indices_[1] = signals_.template readExternalVariableIndex<GensalInputPorts::efd>();
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

    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::setDerivedParams()
    {
      SA_ = 0;
      SB_ = 0;
      if (S12_ != 0)
      {
        RealT s112 = std::sqrt(S10_ / S12_);

        SA_ = (1.2 * s112 + ONE<RealT>) / (s112 + ONE<RealT>);
        SB_ = (1.2 * s112 - ONE<RealT>) / (s112 - ONE<RealT>);
        if (SB_ < SA_)
        {
          SA_ = SB_;
        }
        SB_ = S12_ / ((SA_ - 1.2) * (SA_ - 1.2));
      }
      Xd1_             = Xd_ - Xdp_;
      Xd2_             = Xdp_ - Xl_;
      Xd3_             = (Xdp_ - Xdpp_) / (Xd2_ * Xd2_);
      Xd4_             = (Xdp_ - Xdpp_) / Xd2_;
      Xd5_             = (Xdpp_ - Xl_) / Xd2_;
      Xq2_             = Xq_ - Xdpp_;
      G_               = Ra_ / (Ra_ * Ra_ + Xdpp_ * Xdpp_);
      B_               = -Xdpp_ / (Ra_ * Ra_ + Xdpp_ * Xdpp_);
      va_machine_base_ = mva_base_ * static_cast<RealT>(1.0e6);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
