#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>

#include <GridKit/Model/EMT/Component/Source/Machine/Machine.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/MachineData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase round-rotor synchronous machine
     *
     * System sizes:
     * - Number of equations = 24
     * - Number of independent variables = 24
     */
    template <typename scalar_type, typename index_type>
    Machine<scalar_type, index_type>::Machine()
    {
      size_ = 24;
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Machine<scalar_type, index_type>::Machine(const ModelDataT& data)
      : monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      size_ = 24;
      setDerivedParams();
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Machine<scalar_type, index_type>::~Machine()
    {
    }

    /**
     * @brief Read model parameters from the data object
     */
    template <typename scalar_type, typename index_type>
    void Machine<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      if (data.parameters.contains(Parameter::S))
      {
        S_ = std::get<RealT>(data.parameters.at(Parameter::S));
      }

      if (data.parameters.contains(Parameter::V))
      {
        V_ = std::get<RealT>(data.parameters.at(Parameter::V));
      }

      if (data.parameters.contains(Parameter::f))
      {
        freq_ = std::get<RealT>(data.parameters.at(Parameter::f));
      }

      if (data.parameters.contains(Parameter::H))
      {
        H_ = std::get<RealT>(data.parameters.at(Parameter::H));
      }

      if (data.parameters.contains(Parameter::F))
      {
        Fric_ = std::get<RealT>(data.parameters.at(Parameter::F));
      }

      if (data.parameters.contains(Parameter::Rs))
      {
        Rs_ = std::get<RealT>(data.parameters.at(Parameter::Rs));
      }

      if (data.parameters.contains(Parameter::Ll))
      {
        Ll_ = std::get<RealT>(data.parameters.at(Parameter::Ll));
      }

      if (data.parameters.contains(Parameter::Lmd))
      {
        Lmd_ = std::get<RealT>(data.parameters.at(Parameter::Lmd));
      }

      if (data.parameters.contains(Parameter::Lmq))
      {
        Lmq_ = std::get<RealT>(data.parameters.at(Parameter::Lmq));
      }

      L0_ = Ll_;
      if (data.parameters.contains(Parameter::L0))
      {
        L0_ = std::get<RealT>(data.parameters.at(Parameter::L0));
      }

      if (data.parameters.contains(Parameter::Rfd))
      {
        Rfd_ = std::get<RealT>(data.parameters.at(Parameter::Rfd));
      }

      if (data.parameters.contains(Parameter::Llfd))
      {
        Llfd_ = std::get<RealT>(data.parameters.at(Parameter::Llfd));
      }

      if (data.parameters.contains(Parameter::R1d))
      {
        R1d_ = std::get<RealT>(data.parameters.at(Parameter::R1d));
      }

      if (data.parameters.contains(Parameter::Ll1d))
      {
        Ll1d_ = std::get<RealT>(data.parameters.at(Parameter::Ll1d));
      }

      if (data.parameters.contains(Parameter::R1q))
      {
        R1q_ = std::get<RealT>(data.parameters.at(Parameter::R1q));
      }

      if (data.parameters.contains(Parameter::Ll1q))
      {
        Ll1q_ = std::get<RealT>(data.parameters.at(Parameter::Ll1q));
      }

      if (data.parameters.contains(Parameter::R2q))
      {
        R2q_ = std::get<RealT>(data.parameters.at(Parameter::R2q));
      }

      if (data.parameters.contains(Parameter::Ll2q))
      {
        Ll2q_ = std::get<RealT>(data.parameters.at(Parameter::Ll2q));
      }

      if (data.parameters.contains(Parameter::S10))
      {
        S10_ = std::get<RealT>(data.parameters.at(Parameter::S10));
      }

      if (data.parameters.contains(Parameter::S12))
      {
        S12_ = std::get<RealT>(data.parameters.at(Parameter::S12));
      }

      if (data.parameters.contains(Parameter::p0))
      {
        p0_ = std::get<RealT>(data.parameters.at(Parameter::p0));
      }

      if (data.parameters.contains(Parameter::q0))
      {
        q0_ = std::get<RealT>(data.parameters.at(Parameter::q0));
      }
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename scalar_type, typename index_type>
    void Machine<scalar_type, index_type>::setDerivedParams()
    {
      omega_base_  = TWO<RealT> * std::numbers::pi_v<RealT> * freq_;
      v_peak_base_ = std::sqrt(TWO<RealT> / THREE<RealT>) * V_;
      i_peak_base_ = 0.0;
      if (V_ != 0.0)
      {
        i_peak_base_ = std::numbers::sqrt2_v<RealT> * S_ / (std::sqrt(THREE<RealT>) * V_);
      }
      k_fd_ = 0.0;
      if (Lmd_ != 0.0)
      {
        k_fd_ = Rfd_ / Lmd_;
      }

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
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(static_cast<IdxT>(MachineExternalVariables::MAXIMUM), 3);
      signals_.registerExternalVariableSignals(*this);
      this->setExternalResidualSignal(0, signals_.template getAttachedSignal<MachineExternalVariables::VA>());
      this->setExternalResidualSignal(1, signals_.template getAttachedSignal<MachineExternalVariables::VB>());
      this->setExternalResidualSignal(2, signals_.template getAttachedSignal<MachineExternalVariables::VC>());

      // Set output signals
      if (signals_.template isAssigned<MachineInternalVariables::OMEGA>())
      {
        auto* y = y_.getData();
        signals_.template getSignal<MachineInternalVariables::OMEGA>()->set(&y[1], &(this->getVariableIndex(1)));
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::verify() const
    {
      static constexpr auto EFD = MachineExternalVariables::EFD;
      static constexpr auto PM  = MachineExternalVariables::PM;

      int error_count = 0;

      if (!signals_.template isAttached<MachineExternalVariables::VA>()
          || !signals_.template isAttached<MachineExternalVariables::VB>()
          || !signals_.template isAttached<MachineExternalVariables::VC>())
      {
        Log::error() << "Machine: the bus voltage port is not attached\n";
        ++error_count;
      }

      if (S_ <= 0.0 || V_ <= 0.0 || freq_ <= 0.0 || H_ <= 0.0)
      {
        Log::error() << "Machine: the ratings and inertia must be positive\n";
        ++error_count;
      }

      if (Ll_ <= 0.0 || Lmd_ <= 0.0 || Lmq_ <= 0.0 || L0_ <= 0.0
          || Llfd_ <= 0.0 || Ll1d_ <= 0.0 || Ll1q_ <= 0.0 || Ll2q_ <= 0.0)
      {
        Log::error() << "Machine: the winding inductances must be positive\n";
        ++error_count;
      }

      if (Rfd_ <= 0.0 || R1d_ <= 0.0 || R1q_ <= 0.0 || R2q_ <= 0.0)
      {
        Log::error() << "Machine: the rotor winding resistances must be positive\n";
        ++error_count;
      }

      if (Rs_ < 0.0 || Fric_ < 0.0 || S10_ < 0.0 || S12_ < 0.0)
      {
        Log::error() << "Machine: the stator resistance, friction factor, "
                        "and saturation factors must be nonnegative\n";
        ++error_count;
      }

      if (signals_.template isAttached<EFD>())
      {
        if (!signals_.template isLinked<EFD>())
        {
          Log::error() << "Machine: efd signal attached with no linked exciter\n";
          ++error_count;
        }
      }

      if (signals_.template isAttached<PM>())
      {
        if (!signals_.template isLinked<PM>())
        {
          Log::error() << "Machine: pm signal attached with no linked governor\n";
          ++error_count;
        }
      }

      return error_count;
    }

    /**
     * Initialization of the synchronous machine model
     *
     * The bus terminal instantaneous voltages and the machine power
     * injections are taken as a balanced positive-sequence operating point.
     * All algebra is in machine per unit with peak-value phasors sampled at
     * the initialization instant.
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::initialize()
    {
      using Variables = MachineExternalVariables;

      const RealT pi    = std::numbers::pi_v<RealT>;
      const RealT gamma = TWO<RealT> * pi / THREE<RealT>;

      // Terminal voltage in machine per unit
      const ScalarT va = toMachinePU(signals_.template readExternalVariable<Variables::VA>());
      const ScalarT vb = toMachinePU(signals_.template readExternalVariable<Variables::VB>());
      const ScalarT vc = toMachinePU(signals_.template readExternalVariable<Variables::VC>());

      // Clarke transform gives the rotating peak-value phasor at the
      // initialization instant
      const ScalarT v_re = (TWO<RealT> / THREE<RealT>) *(va - HALF<RealT> * vb - HALF<RealT> * vc);
      const ScalarT v_im = (vb - vc) / std::sqrt(THREE<RealT>);
      const ScalarT vm2  = v_re * v_re + v_im * v_im;

      // Injected current phasor from the scheduled power
      const ScalarT p = static_cast<ScalarT>(p0_) / S_;
      const ScalarT q = static_cast<ScalarT>(q0_) / S_;

      const ScalarT i_re = (p * v_re + q * v_im) / vm2;
      const ScalarT i_im = (p * v_im - q * v_re) / vm2;

      // Saturation from the air-gap flux magnitude behind the leakage
      // impedance
      const ScalarT el_re   = v_re + Rs_ * i_re - Ll_ * i_im;
      const ScalarT el_im   = v_im + Rs_ * i_im + Ll_ * i_re;
      const ScalarT psi_at0 = std::sqrt(el_re * el_re + el_im * el_im);
      const ScalarT ks      = ONE<RealT> / (ONE<RealT> + SB_ * Math::qramp(psi_at0 - SA_));

      const ScalarT Lad = ks * Lmd_;
      const ScalarT Laq = ks * Lmq_;
      const ScalarT Ld  = Lad + Ll_;
      const ScalarT Lq  = Laq + Ll_;

      // Rotor position from the q-axis internal voltage
      const ScalarT eq_re  = v_re + Rs_ * i_re - Lq * i_im;
      const ScalarT eq_im  = v_im + Rs_ * i_im + Lq * i_re;
      const ScalarT theta0 = std::atan2(eq_im, eq_re) - HALF<RealT> * pi;

      const ScalarT ct = std::cos(theta0);
      const ScalarT st = std::sin(theta0);

      // Rotor-frame projections
      const ScalarT id = i_re * ct + i_im * st;
      const ScalarT iq = i_im * ct - i_re * st;
      const ScalarT vd = v_re * ct + v_im * st;
      const ScalarT vq = v_im * ct - v_re * st;

      // Steady stator fluxes with damper currents at zero
      const ScalarT psid = vq + Rs_ * iq;
      const ScalarT psiq = -(vd + Rs_ * id);
      const ScalarT ifd  = (psid + Ld * id) / Lad;

      const ScalarT psifd = (Lad + Llfd_) * ifd - Lad * id;
      const ScalarT psi1d = Lad * (ifd - id);
      const ScalarT psi1q = -Laq * iq;
      const ScalarT psi2q = -Laq * iq;

      const ScalarT te = psid * iq - psiq * id;

      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      y[0]  = theta0;
      y[1]  = ONE<RealT>;
      y[2]  = psid;
      y[3]  = psiq;
      y[4]  = 0.0;
      y[5]  = psifd;
      y[6]  = psi1d;
      y[7]  = psi1q;
      y[8]  = psi2q;
      y[9]  = id;
      y[10] = iq;
      y[11] = 0.0;
      y[12] = ifd;
      y[13] = 0.0;
      y[14] = 0.0;
      y[15] = 0.0;
      y[16] = psid + Ll_ * id;
      y[17] = psiq + Ll_ * iq;
      y[18] = psi_at0;
      y[19] = ks;
      y[20] = te;
      y[21] = id * ct - iq * st;
      y[22] = id * std::cos(theta0 - gamma) - iq * std::sin(theta0 - gamma);
      y[23] = id * std::cos(theta0 + gamma) - iq * std::sin(theta0 + gamma);

      pm_set_ = te + Fric_;
      if (signals_.template isAttached<Variables::PM>())
      {
        signals_.template writeExternalVariable<Variables::PM>(pm_set_);
      }

      efd_set_ = Lmd_ * ifd;
      if (signals_.template isAttached<Variables::EFD>())
      {
        signals_.template writeExternalVariable<Variables::EFD>(efd_set_);
      }

      for (IdxT i = 0; i < size_; ++i)
      {
        yp[static_cast<size_t>(i)] = 0.0;
      }

      // The rotor angle advances at synchronous speed in steady state.
      yp[0] = omega_base_;

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::tagDifferentiable()
    {
      for (size_t j = 0; j < static_cast<size_t>(size_); ++j)
      {
        tag_[j] = j < 9;
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
    int Machine<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Internal residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Machine<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f)
    {
      const RealT pi    = std::numbers::pi_v<RealT>;
      const RealT gamma = TWO<RealT> * pi / THREE<RealT>;

      /* Read variables */
      const ScalarT theta = y[0];
      const ScalarT omega = y[1];
      const ScalarT psid  = y[2];
      const ScalarT psiq  = y[3];
      const ScalarT psi0  = y[4];
      const ScalarT psifd = y[5];
      const ScalarT psi1d = y[6];
      const ScalarT psi1q = y[7];
      const ScalarT psi2q = y[8];
      const ScalarT id    = y[9];
      const ScalarT iq    = y[10];
      const ScalarT i0    = y[11];
      const ScalarT ifd   = y[12];
      const ScalarT i1d   = y[13];
      const ScalarT i1q   = y[14];
      const ScalarT i2q   = y[15];
      const ScalarT psiad = y[16];
      const ScalarT psiaq = y[17];
      const ScalarT psiat = y[18];
      const ScalarT ks    = y[19];
      const ScalarT te    = y[20];
      const ScalarT isa   = y[21];
      const ScalarT isb   = y[22];
      const ScalarT isc   = y[23];

      /* Read derivatives */
      const ScalarT theta_dot = yp[0];
      const ScalarT omega_dot = yp[1];
      const ScalarT psid_dot  = yp[2];
      const ScalarT psiq_dot  = yp[3];
      const ScalarT psi0_dot  = yp[4];
      const ScalarT psifd_dot = yp[5];
      const ScalarT psi1d_dot = yp[6];
      const ScalarT psi1q_dot = yp[7];
      const ScalarT psi2q_dot = yp[8];

      // Set coupling variable aliases
      const ScalarT va_pu = y_ext[0] / v_peak_base_;
      const ScalarT vb_pu = y_ext[1] / v_peak_base_;
      const ScalarT vc_pu = y_ext[2] / v_peak_base_;
      const ScalarT efd   = k_fd_ * y_ext[3];
      const ScalarT pm    = y_ext[4];

      const ScalarT ct_a = std::cos(theta);
      const ScalarT st_a = std::sin(theta);
      const ScalarT ct_b = std::cos(theta - gamma);
      const ScalarT st_b = std::sin(theta - gamma);
      const ScalarT ct_c = std::cos(theta + gamma);
      const ScalarT st_c = std::sin(theta + gamma);

      const ScalarT vd = (TWO<RealT> / THREE<RealT>) *(va_pu * ct_a + vb_pu * ct_b + vc_pu * ct_c);
      const ScalarT vq = -(TWO<RealT> / THREE<RealT>) *(va_pu * st_a + vb_pu * st_b + vc_pu * st_c);
      const ScalarT v0 = (va_pu + vb_pu + vc_pu) / THREE<RealT>;

      const ScalarT Lad = ks * Lmd_;
      const ScalarT Laq = ks * Lmq_;

      const RealT inv_omega_base = ONE<RealT> / omega_base_;

      /* 9 machine differential equations */
      f[0] = theta_dot - omega_base_ * omega;
      f[1] = TWO<RealT> * H_ * omega_dot - pm / omega + te + Fric_ * omega;
      f[2] = inv_omega_base * psid_dot - omega * psiq - Rs_ * id - vd;
      f[3] = inv_omega_base * psiq_dot + omega * psid - Rs_ * iq - vq;
      f[4] = inv_omega_base * psi0_dot - Rs_ * i0 - v0;
      f[5] = inv_omega_base * psifd_dot + Rfd_ * ifd - efd;
      f[6] = inv_omega_base * psi1d_dot + R1d_ * i1d;
      f[7] = inv_omega_base * psi1q_dot + R1q_ * i1q;
      f[8] = inv_omega_base * psi2q_dot + R2q_ * i2q;

      /* 15 machine algebraic equations */
      f[9]  = psid + (Lad + Ll_) * id - Lad * (ifd + i1d);
      f[10] = psiq + (Laq + Ll_) * iq - Laq * (i1q + i2q);
      f[11] = psi0 + L0_ * i0;
      f[12] = psifd - (Lad + Llfd_) * ifd - Lad * i1d + Lad * id;
      f[13] = psi1d - Lad * ifd - (Lad + Ll1d_) * i1d + Lad * id;
      f[14] = psi1q - (Laq + Ll1q_) * i1q - Laq * i2q + Laq * iq;
      f[15] = psi2q - Laq * i1q - (Laq + Ll2q_) * i2q + Laq * iq;
      f[16] = psiad - psid - Ll_ * id;
      f[17] = psiaq - psiq - Ll_ * iq;
      f[18] = psiat - std::sqrt(psiad * psiad + psiaq * psiaq);
      f[19] = ks * (ONE<RealT> + SB_ * Math::qramp(psiat - SA_)) - ONE<RealT>;
      f[20] = te - (psid * iq - psiq * id);
      f[21] = isa - (id * ct_a - iq * st_a + i0);
      f[22] = isb - (id * ct_b - iq * st_b + i0);
      f[23] = isc - (id * ct_c - iq * st_c + i0);

      return 0;
    }

    /**
     * @brief External residual
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Machine<scalar_type, index_type>::evaluateExternalResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        [[maybe_unused]] const ScalarT* y_ext,
        [[maybe_unused]] const ScalarT* yp_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT isa = y[21];
      const ScalarT isb = y[22];
      const ScalarT isc = y[23];

      f_ext[0] = toSystemSI(isa);
      f_ext[1] = toSystemSI(isb);
      f_ext[2] = toSystemSI(isc);

      return 0;
    }

    /**
     * @brief Gather external variables and index maps.
     *
     * The latched setpoints back the field voltage and mechanical power
     * inputs when no controller is attached; the base gather then refreshes
     * every attached slot through its signal.
     */
    template <typename scalar_type, typename index_type>
    void Machine<scalar_type, index_type>::gatherExternalVariables()
    {
      y_ext_[static_cast<size_t>(MachineExternalVariables::EFD)] = efd_set_;
      y_ext_[static_cast<size_t>(MachineExternalVariables::PM)]  = pm_set_;

      Component<scalar_type, index_type>::gatherExternalVariables();
    }

    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      evaluateInternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f);
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief External residual contributions to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::evaluateExternalResidual()
    {
      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();

      return 0;
    }

    /**
     * @brief Residual contribution of the machine is pushed to the bus.
     *
     */
    template <typename scalar_type, typename index_type>
    int Machine<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Machine<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Machine<scalar_type, index_type>::initializeMonitor()
    {
      using Variable  = typename ModelDataT::MonitorableVariables;
      using Variables = MachineExternalVariables;

      monitor_->set(Variable::theta, [this]
                    { return y_.getData()[0]; });
      monitor_->set(Variable::omega, [this]
                    { return y_.getData()[1]; });
      monitor_->set(Variable::te, [this]
                    { return y_.getData()[20]; });
      monitor_->set(Variable::ifd, [this]
                    { return y_.getData()[12]; });
      monitor_->set(Variable::efd, [this]
                    { return signals_.template isAttached<MachineExternalVariables::EFD>()
                                 ? signals_.template readExternalVariable<MachineExternalVariables::EFD>()
                                 : efd_set_; });
      monitor_->set(Variable::ks, [this]
                    { return y_.getData()[19]; });
      monitor_->set(Variable::psi_at, [this]
                    { return y_.getData()[18]; });
      monitor_->set(Variable::ia, [this]
                    { return toSystemSI(y_.getData()[21]); });
      monitor_->set(Variable::ib, [this]
                    { return toSystemSI(y_.getData()[22]); });
      monitor_->set(Variable::ic, [this]
                    { return toSystemSI(y_.getData()[23]); });
      monitor_->set(Variable::p, [this]
                    {
                      const auto* y = y_.getData();
                      return signals_.template readExternalVariable<Variables::VA>() * toSystemSI(y[21])
                             + signals_.template readExternalVariable<Variables::VB>() * toSystemSI(y[22])
                             + signals_.template readExternalVariable<Variables::VC>() * toSystemSI(y[23]); });
      monitor_->set(Variable::q, [this]
                    {
                      // Reactive power from the rotor-frame quantities scaled
                      // to the machine base
                      const auto* y  = y_.getData();
                      const auto  ct = std::cos(y[0]);
                      const auto  st = std::sin(y[0]);
                      const auto  va = toMachinePU(signals_.template readExternalVariable<Variables::VA>());
                      const auto  vb = toMachinePU(signals_.template readExternalVariable<Variables::VB>());
                      const auto  vc = toMachinePU(signals_.template readExternalVariable<Variables::VC>());
                      const auto  gamma = TWO<RealT> * std::numbers::pi_v<RealT> / THREE<RealT>;
                      const auto  vd = (TWO<RealT> / THREE<RealT>)
                                      * (va * ct + vb * std::cos(y[0] - gamma) + vc * std::cos(y[0] + gamma));
                      const auto vq = -(TWO<RealT> / THREE<RealT>)
                                      * (va * st + vb * std::sin(y[0] - gamma) + vc * std::sin(y[0] + gamma));
                      return S_ * (vq * y[9] - vd * y[10]); });
      monitor_->set(Variable::id, [this]
                    { return y_.getData()[9]; });
      monitor_->set(Variable::iq, [this]
                    { return y_.getData()[10]; });
    }

  } // namespace EMT
} // namespace GridKit
