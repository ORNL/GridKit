/**
 * @file ReecbImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REECB electrical-control model.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /// Logger used for REECB diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct REECB with its documented parameter defaults
       *
       * The terminal bus is retained, the model is sized, and no monitor or
       * signal connection is created.
       *
       * @param[in] bus Terminal bus measured by the controller.
       */
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      /**
       * @brief Construct REECB from model data
       *
       * @param[in] bus Terminal bus measured by the controller.
       * @param[in] data Model parameters and monitor selections.
       */
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(ReecbInternalVariables::MAXIMUM);
      }

      /**
       * @brief Destroy the electrical controller and its optional variable monitor.
       */
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::~Reecb()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate model vectors and wire assigned current-command outputs
       *
       * Sizes the state, residual, bus, and signal-interface buffers, initializes
       * identity index maps, and points assigned command nodes at the internal
       * system-base states that REECB publishes. Repeated allocation reuses the
       * existing model vectors and signal links.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::allocate()
      {
        const auto IQCMD = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        const auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        const auto signal_size = static_cast<size_t>(ReecbExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (signals_.template isAssigned<ReecbInternalVariables::IQCMD>())
        {
          signals_.template getSignalNode<ReecbInternalVariables::IQCMD>()->set(
              &y[IQCMD],
              &(this->getVariableIndex(static_cast<IdxT>(IQCMD))));
        }

        if (signals_.template isAssigned<ReecbInternalVariables::IPCMD>())
        {
          signals_.template getSignalNode<ReecbInternalVariables::IPCMD>()->set(
              &y[IPCMD],
              &(this->getVariableIndex(static_cast<IdxT>(IPCMD))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the REECB configuration
       *
       * Checks parameter-loading errors, finiteness and static relationships,
       * system/component bases and conversion ratios, the terminal bus, and
       * attached optional signals. Operating-point feasibility is checked by
       * initialize().
       *
       * @return Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        checkConfiguration(bus_ != nullptr, "terminal bus is required", ret);

        const RealT component_power_base = componentPowerBase();
        const bool  valid_component_base = std::isfinite(component_power_base) && component_power_base > ZERO<RealT>;
        const bool  valid_system_base    = std::isfinite(va_system_base_) && va_system_base_ > ZERO<RealT>;
        checkConfiguration(valid_component_base, "component power base must be finite and positive", ret);
        checkConfiguration(valid_system_base, "system power base must be finite and positive", ret);
        if (valid_component_base && valid_system_base)
        {
          const RealT system_to_component = va_system_base_ / component_power_base;
          const RealT component_to_system = component_power_base / va_system_base_;
          checkConfiguration(
              std::isfinite(system_to_component)
                  && system_to_component > ZERO<RealT>
                  && std::isfinite(component_to_system)
                  && component_to_system > ZERO<RealT>,
              "system/component power-base conversion ratios must be finite and positive",
              ret);
        }

        checkConfiguration(std::isfinite(Trv_), "Trv must be finite", ret);
        checkConfiguration(std::isfinite(Tp_), "Tp must be finite", ret);
        checkConfiguration(std::isfinite(Vref0_), "Vref0 must be finite", ret);

        const bool finite_voltage_thresholds = std::isfinite(Vdip_) && std::isfinite(Vup_);
        checkConfiguration(finite_voltage_thresholds, "Vdip and Vup must be finite", ret);
        if (finite_voltage_thresholds)
        {
          checkConfiguration(Vdip_ < Vup_, "Vdip must be less than Vup", ret);
        }

        const bool finite_voltage_deadband = std::isfinite(dbd1_) && std::isfinite(dbd2_);
        checkConfiguration(finite_voltage_deadband, "dbd1 and dbd2 must be finite", ret);
        if (finite_voltage_deadband)
        {
          checkConfiguration(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required", ret);
        }

        checkConfiguration(std::isfinite(kqv_), "kqv must be finite", ret);

        const bool finite_injection_limits = std::isfinite(Iql1_) && std::isfinite(Iqh1_);
        checkConfiguration(finite_injection_limits, "Iql1 and Iqh1 must be finite", ret);
        if (finite_injection_limits)
        {
          checkConfiguration(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1", ret);
        }

        const bool finite_reactive_limits = std::isfinite(Qmin_) && std::isfinite(Qmax_);
        checkConfiguration(finite_reactive_limits, "Qmin and Qmax must be finite", ret);
        if (finite_reactive_limits)
        {
          checkConfiguration(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax", ret);
        }

        checkConfiguration(std::isfinite(Kqp_), "Kqp must be finite", ret);
        checkConfiguration(std::isfinite(Kqi_), "Kqi must be finite", ret);

        const bool finite_voltage_limits = std::isfinite(Vmin_) && std::isfinite(Vmax_);
        checkConfiguration(finite_voltage_limits, "Vmin and Vmax must be finite", ret);
        if (finite_voltage_limits)
        {
          checkConfiguration(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax", ret);
        }

        checkConfiguration(std::isfinite(Kvp_), "Kvp must be finite", ret);
        checkConfiguration(std::isfinite(Kvi_), "Kvi must be finite", ret);
        checkConfiguration(std::isfinite(Tiq_), "Tiq must be finite", ret);
        checkConfiguration(std::isfinite(Tpord_), "Tpord must be finite", ret);

        const bool finite_ramp_limits = std::isfinite(dPmin_) && std::isfinite(dPmax_);
        checkConfiguration(finite_ramp_limits, "dPmin and dPmax must be finite", ret);
        if (finite_ramp_limits)
        {
          checkConfiguration(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required", ret);
        }

        const bool finite_active_limits = std::isfinite(Pmin_) && std::isfinite(Pmax_);
        checkConfiguration(finite_active_limits, "Pmin and Pmax must be finite", ret);
        if (finite_active_limits)
        {
          checkConfiguration(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax", ret);
        }

        checkConfiguration(std::isfinite(Imax_) && Imax_ > ZERO<RealT>, "Imax must be finite and positive", ret);

        checkConfiguration(!(PfFlag_ && QFlag_ && !VFlag_),
                           "power-factor control cannot drive the direct-voltage reference (PfFlag = 1 with QFlag = 1, VFlag = 0)",
                           ret);

        auto check_optional_signal = [&]<ReecbExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>() && !signals_.template isLinked<variable>())
          {
            Log::error() << "Reecb: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_optional_signal.template operator()<ReecbExternalVariables::PE>("pe");
        check_optional_signal.template operator()<ReecbExternalVariables::QGEN>("qgen");
        check_optional_signal.template operator()<ReecbExternalVariables::QEXT>("qext");
        check_optional_signal.template operator()<ReecbExternalVariables::PFAREF>("pfaref");
        check_optional_signal.template operator()<ReecbExternalVariables::PREF>("pref");

        return ret;
      }

      /**
       * @brief Initialize REECB from the initial current commands and feedback
       *
       * Preserves the system-base command states, consumes attached initialized
       * active/reactive-power feedback or reconstructs unattached feedback, and
       * constructs the remaining states and reference setpoints by forward
       * evaluation of the residual expressions, so an admissible operating
       * point starts at a steady state.
       *
       * @pre allocate() has completed.
       * @pre verify() reports a valid parameter and port configuration.
       * @pre The terminal bus and current-command states are initialized.
       *
       * @post On failure no state, derivative, parameter, or signal storage is
       *       modified.
       *
       * @return 0 on success; nonzero when an allocation, configuration,
       *         initial-value, current-circle, or limiter-interior check fails.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initialize()
      {
        const auto VMEAS = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ  = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV  = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV    = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD  = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT    = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto ILMAX = static_cast<size_t>(ReecbInternalVariables::ILMAX);
        const auto IQCMD = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        if (!allocated_)
        {
          Log::error() << "Reecb: allocate must complete before initialize\n";
          return 1;
        }

        if (verify() > 0)
        {
          Log::error() << "Reecb: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();

        const RealT ipcmd0_system = static_cast<RealT>(y[IPCMD]);
        const RealT iqcmd0_system = static_cast<RealT>(y[IQCMD]);
        const RealT ipcmd0        = toComponentBase(ipcmd0_system);
        const RealT iqcmd0        = toComponentBase(iqcmd0_system);
        const RealT vr0           = static_cast<RealT>(Vr());
        const RealT vi0           = static_cast<RealT>(Vi());
        const RealT vt0           = std::sqrt(vr0 * vr0 + vi0 * vi0);
        const RealT vmeas0        = vt0;
        const RealT vmeas_safe0   = Math::max(vmeas0, VMEAS_MINIMUM);

        RealT pe0_system   = toSystemBase(ipcmd0 * vmeas_safe0);
        RealT qgen0_system = toSystemBase(iqcmd0 * vmeas_safe0);

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          pe0_system = static_cast<RealT>(
              signals_.template readExternalVariable<ReecbExternalVariables::PE>());
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          qgen0_system = static_cast<RealT>(
              signals_.template readExternalVariable<ReecbExternalVariables::QGEN>());
        }

        const RealT pmeas0 = toComponentBase(pe0_system);
        const RealT qgen0  = toComponentBase(qgen0_system);
        RealT       vref0  = vmeas0;
        if (Vref0_given_)
        {
          vref0 = Vref0_;
        }

        if (!std::isfinite(vr0) || !std::isfinite(vi0) || !std::isfinite(vt0) || !std::isfinite(vmeas_safe0)
            || !std::isfinite(ipcmd0) || !std::isfinite(iqcmd0) || !std::isfinite(pmeas0)
            || !std::isfinite(qgen0) || !std::isfinite(vref0))
        {
          Log::error() << "Reecb: initial bus, command, and feedback values must be finite\n";
          return 1;
        }
        if (vt0 <= ZERO<RealT>)
        {
          Log::error() << "Reecb: initial terminal-voltage magnitude must be positive\n";
          return 1;
        }

        // Mirrors of the residual limiter chain at the initial point.
        const RealT verr0         = Math::deadband2(vref0 - vmeas0, dbd1_, dbd2_);
        const RealT iqv0          = Math::clamp(kqv_ * verr0, Iql1_, Iqh1_);
        const RealT ilmax_squared = Imax_ * Imax_ - pq_on_ * ipcmd0 * ipcmd0 - pq_off_ * iqcmd0 * iqcmd0;

        if (!std::isfinite(ilmax_squared) || ilmax_squared <= ZERO<RealT>)
        {
          Log::error() << "Reecb: initial operating point leaves no low-priority current capacity\n";
          return 1;
        }

        const RealT ilmax0 = std::sqrt(ilmax_squared);
        const RealT iqmax0 = pq_on_ * ilmax0 + pq_off_ * Imax_;
        const RealT ipmax0 = pq_on_ * Imax_ + pq_off_ * ilmax0;

        if (ipcmd0 <= ZERO<RealT> || ipcmd0 >= ipmax0 || iqcmd0 <= -iqmax0 || iqcmd0 >= iqmax0)
        {
          Log::error() << "Reecb: initial current commands must lie strictly inside their limiter ranges\n";
          return 1;
        }

        // The algebraic command rows reproduce their limiter outputs through
        // the smooth-clamp inverse. A command no input can produce leaves the
        // inverse nonfinite, which the finiteness test below rejects.
        const RealT ipraw0 = unclamp(ipcmd0, ZERO<RealT>, ipmax0);
        const RealT iqraw0 = unclamp(iqcmd0, -iqmax0, iqmax0);
        const RealT iqctl0 = iqraw0 - iqv0;
        const RealT pord0  = vmeas_safe0 * ipraw0;

        const RealT pref0_system = toSystemBase(pord0);

        if (pord0 < Pmin_ || pord0 > Pmax_)
        {
          Log::error() << "Reecb: recovered active-power order is outside Pmin/Pmax\n";
          return 1;
        }

        // An integrating path holds its feedback only where the clamp can
        // reproduce it: strictly inside the limits, or collapsed onto it.
        auto reproducible = [](RealT value, RealT lower, RealT upper)
        {
          return (lower < value && value < upper) || (lower == upper && value == lower);
        };

        if (q_pi_on_ * Kqi_ != ZERO<RealT> && !reproducible(qgen0, Qmin_, Qmax_))
        {
          Log::error() << "Reecb: reactive-power integral path is not at equilibrium\n";
          return 1;
        }
        if (q_pi_on_ * Kvi_ != ZERO<RealT> && !reproducible(vmeas0, Vmin_, Vmax_))
        {
          Log::error() << "Reecb: voltage-control integral path is not at equilibrium\n";
          return 1;
        }

        // The reactive channel reproduces the power feedback when the reactive
        // PI is enabled, and the reactive-current command otherwise. Collapsed
        // limits already pin the clamp output onto the feedback.
        RealT qtarget0 = ZERO<RealT>;
        if (!QFlag_)
        {
          qtarget0 = iqctl0 * vmeas_safe0;
        }
        else if (VFlag_ && Qmin_ < qgen0 && qgen0 < Qmax_)
        {
          qtarget0 = unclamp(qgen0, Qmin_, Qmax_);
        }

        RealT qref0      = ZERO<RealT>;
        RealT qext0_port = ZERO<RealT>;
        RealT pfaref0    = ZERO<RealT>;

        if (QFlag_ && !VFlag_)
        {
          // Direct-voltage mode publishes the measurement as its reference.
          qext0_port = vmeas0;
        }
        else if (PfFlag_)
        {
          if (pmeas0 == ZERO<RealT> && qtarget0 != ZERO<RealT>)
          {
            Log::error() << "Reecb: power-factor mode cannot reproduce the reactive target at zero active power\n";
            return 1;
          }
          if (pmeas0 != ZERO<RealT>)
          {
            pfaref0 = std::atan(qtarget0 / pmeas0);
          }
          qref0 = pmeas0 * std::tan(pfaref0);

          // Angle resolution collapses toward the tangent pole, so the
          // published angle must still carry its own target back.
          if (std::abs(qref0 - qtarget0) > std::abs(qtarget0) * INITIALIZATION_TOLERANCE)
          {
            Log::error() << "Reecb: power-factor angle cannot reproduce the reactive target\n";
            return 1;
          }
          qext0_port = toSystemBase(qref0);
        }
        else
        {
          qext0_port = toSystemBase(qtarget0);
          qref0      = toComponentBase(qext0_port);
        }

        const RealT eq0 = Math::clamp(qref0, Qmin_, Qmax_) - qgen0;

        // The Q-PI order carries the measurement the voltage channel
        // subtracts; collapsed limits pin the clamp output there already.
        RealT xpiq0 = ZERO<RealT>;
        if (QFlag_ && VFlag_ && Vmin_ < vmeas0 && vmeas0 < Vmax_)
        {
          xpiq0 = unclamp(vmeas0, Vmin_, Vmax_) - Kqp_ * eq0;
        }

        const RealT vpiq0 = Math::clamp(Kqp_ * eq0 + xpiq0, Vmin_, Vmax_);
        const RealT epiv0 = q_pi_on_ * vpiq0 + v_ref_on_ * qext0_port - q_on_ * vmeas0;

        RealT qv0   = ZERO<RealT>;
        RealT xpiv0 = ZERO<RealT>;

        if (QFlag_)
        {
          if (iqctl0 <= -iqmax0 || iqctl0 >= iqmax0)
          {
            Log::error() << "Reecb: initial voltage-controller current is outside its limiter range\n";
            return 1;
          }
          xpiv0 = unclamp(iqctl0, -iqmax0, iqmax0) - Kvp_ * epiv0;
        }
        else
        {
          // The lag state carries the same quotient the QV row forms.
          qv0 = qref0 / vmeas_safe0;
        }

        if (!std::isfinite(verr0) || !std::isfinite(iqv0) || !std::isfinite(ilmax0)
            || !std::isfinite(iqmax0) || !std::isfinite(ipmax0) || !std::isfinite(ipraw0)
            || !std::isfinite(iqraw0) || !std::isfinite(pord0) || !std::isfinite(pref0_system)
            || !std::isfinite(qref0) || !std::isfinite(qext0_port) || !std::isfinite(pfaref0)
            || !std::isfinite(eq0) || !std::isfinite(xpiq0) || !std::isfinite(epiv0)
            || !std::isfinite(xpiv0) || !std::isfinite(qv0))
        {
          Log::error() << "Reecb: initialization produced a nonfinite value\n";
          return 1;
        }

        y[VMEAS] = vmeas0;
        y[PMEAS] = pmeas0;
        y[XPIQ]  = xpiq0;
        y[XPIV]  = xpiv0;
        y[QV]    = qv0;
        y[PORD]  = pord0;
        y[VT]    = vt0;
        y[ILMAX] = ilmax0;

        if (!Vref0_given_)
        {
          Vref0_ = vref0;
        }

        pe_set_     = static_cast<ScalarT>(pe0_system);
        qgen_set_   = static_cast<ScalarT>(qgen0_system);
        qext_set_   = static_cast<ScalarT>(qext0_port);
        pfaref_set_ = static_cast<ScalarT>(pfaref0);
        pref_set_   = static_cast<ScalarT>(pref0_system);

        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::QEXT>(qext_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PFAREF>(pfaref_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PREF>(pref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The two measurement filters, two PI states, reactive-current lag, and
       * active-power order carry derivatives; all other rows are algebraic.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::tagDifferentiable()
      {
        const auto VMEAS = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ  = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV  = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV    = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD  = static_cast<size_t>(ReecbInternalVariables::PORD);

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[VMEAS] = true;
        tag_[PMEAS] = true;
        tag_[XPIQ]  = true;
        tag_[XPIV]  = true;
        tag_[QV]    = true;
        tag_[PORD]  = true;
        return 0;
      }

      /**
       * @brief Set one absolute tolerance for every REECB variable
       *
       * @param[in] rel_tol Solver relative tolerance used as the absolute floor.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Evaluate the model residuals
       *
       * Starts from latched values, refreshes attached signals and their indices,
       * refreshes terminal-bus voltage, and evaluates the internal residual.
       * REECB contributes no bus residual.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateResidual()
      {
        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        ws_[PE]     = pe_set_;
        ws_[QGEN]   = qgen_set_;
        ws_[QEXT]   = qext_set_;
        ws_[PFAREF] = pfaref_set_;
        ws_[PREF]   = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          ws_[PE] = signals_.template readExternalVariable<ReecbExternalVariables::PE>();
          ws_indices_[PE] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PE>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          ws_[QGEN] = signals_.template readExternalVariable<ReecbExternalVariables::QGEN>();
          ws_indices_[QGEN] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::QGEN>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          ws_[QEXT] = signals_.template readExternalVariable<ReecbExternalVariables::QEXT>();
          ws_indices_[QEXT] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::QEXT>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          ws_[PFAREF] =
              signals_.template readExternalVariable<ReecbExternalVariables::PFAREF>();
          ws_indices_[PFAREF] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PFAREF>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<ReecbExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PREF>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        evaluateInternalResidual(y_.getData(), yp_.getData(), wb_.data(), ws_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the REECB signal interface
       *
       * @return Interface used to assign current-command outputs and attach
       *         optional feedback and reference signals.
       */
      template <typename scalar_type, typename index_type>
      auto Reecb<scalar_type, index_type>::getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              ReecbInternalVariables,
                              ReecbExternalVariables>&
      {
        return signals_;
      }

      /**
       * @brief Access the optional variable monitor
       *
       * @return Monitor, or nullptr when constructed without model data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Reecb<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the REECB internal residual
       *
       * The branch-free equation body preserves a fixed dependency structure;
       * parameter-selected paths enter through selector masks resolved by
       * setDerivedParameters(). The `ILMAX` row uses a signed-square
       * continuation, while its magnitude supplies the limiter bounds, so a
       * negative nonlinear iterate does not invert either range.
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] wb Terminal-bus voltage components.
       * @param[in] ws External signal values in their documented port units and bases.
       * @param[out] f Internal residuals.
       * @pre Jacobian evaluation requires `y[ILMAX] != 0`; initialization
       *      rejects the zero-capacity point.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline int
      Reecb<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto VMEAS = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ  = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV  = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV    = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD  = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT    = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto ILMAX = static_cast<size_t>(ReecbInternalVariables::ILMAX);
        const auto IQCMD = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        const ScalarT vmeas        = y[VMEAS];
        const ScalarT pmeas        = y[PMEAS];
        const ScalarT xpiq         = y[XPIQ];
        const ScalarT xpiv         = y[XPIV];
        const ScalarT qv           = y[QV];
        const ScalarT pord         = y[PORD];
        const ScalarT vt           = y[VT];
        const ScalarT ilmax        = y[ILMAX];
        const ScalarT iqcmd_system = y[IQCMD];
        const ScalarT ipcmd_system = y[IPCMD];

        const ScalarT vmeas_dot = yp[VMEAS];
        const ScalarT pmeas_dot = yp[PMEAS];
        const ScalarT xpiq_dot  = yp[XPIQ];
        const ScalarT xpiv_dot  = yp[XPIV];
        const ScalarT qv_dot    = yp[QV];
        const ScalarT pord_dot  = yp[PORD];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT pe     = toComponentBase(ws[PE]);
        const ScalarT qgen   = toComponentBase(ws[QGEN]);
        const ScalarT extref = ws[QEXT];
        const ScalarT pfaref = ws[PFAREF];
        const ScalarT pref   = toComponentBase(ws[PREF]);
        const ScalarT iqcmd  = toComponentBase(iqcmd_system);
        const ScalarT ipcmd  = toComponentBase(ipcmd_system);

        const ScalarT vmeas_safe = Math::max(vmeas, VMEAS_MINIMUM);
        const ScalarT sdip       = Math::inside(vt, Vdip_, Vup_);
        const ScalarT verr       = Math::deadband2(Vref0_ - vmeas, dbd1_, dbd2_);
        const ScalarT iqv        = Math::clamp(kqv_ * verr, Iql1_, Iqh1_);
        // The Volt/VAr channel is a system-base reactive power unless
        // direct-voltage mode selects it as a terminal-voltage reference,
        // which takes no power-base conversion.
        const ScalarT qref       = q_ref_on_ * (pf_on_ * pmeas * std::tan(pfaref) + pf_off_ * toComponentBase(extref));
        const ScalarT eq         = Math::clamp(qref, Qmin_, Qmax_) - qgen;
        const ScalarT vpiq       = Math::clamp(Kqp_ * eq + xpiq, Vmin_, Vmax_);
        const ScalarT epiv       = q_pi_on_ * vpiq + v_ref_on_ * extref - q_on_ * vmeas;
        const ScalarT fpord      = (pref - pord) / Tpord_;
        const ScalarT rpord      = aslew(fpord, dPmin_, dPmax_);
        const ScalarT ilcap      = std::sqrt(ilmax * ilmax);
        const ScalarT iqmax      = pq_on_ * ilcap + pq_off_ * Imax_;
        const ScalarT ipmax      = pq_on_ * Imax_ + pq_off_ * ilcap;
        const ScalarT iqbase     = Math::clamp(Kvp_ * epiv + xpiv, -iqmax, iqmax);
        const ScalarT iqraw      = q_on_ * iqbase + q_off_ * qv + iqv;

        f[VMEAS] = -vmeas_dot + (vt - vmeas) / Trv_;
        f[PMEAS] = -pmeas_dot + (pe - pmeas) / Tp_;
        f[XPIQ]  = -xpiq_dot + q_pi_on_ * sdip * Math::antiwindup(Kqp_ * eq + xpiq, Kqi_ * eq, Vmin_, Vmax_);
        f[XPIV]  = -xpiv_dot + q_on_ * sdip * awband(Kvp_ * epiv + xpiv, Kvi_ * epiv, iqmax);
        f[QV]    = -qv_dot + q_off_ * sdip * (qref / vmeas_safe - qv) / Tiq_;
        f[PORD]  = -pord_dot + sdip * Math::antiwindup(pord, rpord, Pmin_, Pmax_);
        f[VT]    = -vt * vt + vr * vr + vi * vi;
        f[ILMAX] = -ilmax * ilcap + Imax_ * Imax_ - pq_on_ * ipcmd * ipcmd - pq_off_ * iqcmd * iqcmd;
        f[IQCMD] = -iqcmd + Math::clamp(iqraw, -iqmax, iqmax);
        f[IPCMD] = -ipcmd + Math::clamp(pord / vmeas_safe, ZERO<RealT>, ipmax);

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Smooth asymmetric slew-rate limiter
       *
       * @param[in] f Unconstrained rate.
       * @param[in] rate_min Negative rate limit.
       * @param[in] rate_max Positive rate limit.
       * @return Limited rate.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Reecb<scalar_type, index_type>::aslew(
          const ScalarT f,
          const RealT   rate_min,
          const RealT   rate_max)
      {
        assert(rate_min < ZERO<RealT> && ZERO<RealT> < rate_max);

        return f
               / (ONE<RealT>
                  + Math::ramp(f / rate_max - ONE<RealT>)
                  + Math::ramp(f / rate_min - ONE<RealT>));
      }

      /**
       * @brief Smooth anti-windup derivative within a moving symmetric band
       *
       * Math::antiwindup over [-band, band] with a band edge that is an
       * algebraic quantity, so differentiation carries the band's own
       * contributions through the gate.
       *
       * @param[in] x Limited PI state.
       * @param[in] f Pre-limit derivative of x.
       * @param[in] band Nonnegative symmetric band edge.
       * @return Anti-windup-limited derivative.
       *
       * @todo Fold moving-limit support into Math::antiwindup in CommonMath.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Reecb<scalar_type, index_type>::awband(
          const ScalarT x,
          const ScalarT f,
          const ScalarT band)
      {
        const ScalarT above_min = Math::sigmoid(x + band);
        const ScalarT below_max = Math::sigmoid(band - x);

        return (above_min * below_max                          //
                + (ONE<RealT> - below_max) * Math::sigmoid(-f) //
                + (ONE<RealT> - above_min) * Math::sigmoid(f))
               * f;
      }

      /**
       * @brief Record one failed configuration condition
       *
       * @param[in] condition Required condition.
       * @param[in] message Error message when `condition` is false.
       * @param[in,out] errors Accumulated configuration-error count.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::checkConfiguration(
          bool condition, const char* message, int& errors)
      {
        if (!condition)
        {
          Log::error() << "Reecb: " << message << '\n';
          errors += 1;
        }
      }

      /**
       * @brief Load one real-valued parameter
       *
       * Real and integer serialized values are accepted. Any other stored type
       * records a loading error while preserving the existing value.
       *
       * @param[in] data Model parameter data.
       * @param[in] parameter Parameter key to load.
       * @param[in,out] target Stored parameter value.
       * @param[in] name Serialized parameter name for diagnostics.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::loadRealParameter(
          const ModelDataT& data,
          ReecbParameters   parameter,
          RealT&            target,
          const char*       name)
      {
        if (!data.parameters.contains(parameter))
        {
          return;
        }

        const auto& value = data.parameters.at(parameter);
        if (const auto* real_value = std::get_if<RealT>(&value))
        {
          target = *real_value;
        }
        else if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          target = static_cast<RealT>(*index_value);
        }
        else
        {
          Log::error() << "Reecb: parameter '" << name << "' must be numeric\n";
          ++parameter_error_count_;
        }
      }

      /**
       * @brief Load one optional Boolean parameter
       *
       * Any non-Boolean stored type records a loading error while preserving
       * the existing default.
       *
       * @param[in] data Model parameter data.
       * @param[in] parameter Parameter key to load.
       * @param[in,out] target Stored Boolean value.
       * @param[in] name Serialized parameter name for diagnostics.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::loadBooleanParameter(
          const ModelDataT& data,
          ReecbParameters   parameter,
          bool&             target,
          const char*       name)
      {
        if (!data.parameters.contains(parameter))
        {
          return;
        }

        const auto& value = data.parameters.at(parameter);
        if (const auto* bool_value = std::get_if<bool>(&value))
        {
          target = *bool_value;
        }
        else
        {
          Log::error() << "Reecb: parameter '" << name << "' must be boolean\n";
          ++parameter_error_count_;
        }
      }

      /**
       * @brief Validate and floor one explicit controller lag
       *
       * Nonfinite and negative values record errors before replacement so
       * verify() retains the evidence. A valid value below 1 ms is raised and
       * reported through the return value, preserving an explicit Hessenberg
       * residual and avoiding division by zero.
       *
       * @param[in,out] value Time constant to validate and floor.
       * @param[in] name Parameter name for diagnostics.
       * @return true only when a valid nonnegative value was raised.
       */
      template <typename scalar_type, typename index_type>
      bool Reecb<scalar_type, index_type>::floorTimeConstant(
          RealT& value, const char* name)
      {
        if (!std::isfinite(value))
        {
          Log::error() << "Reecb: " << name << " must be finite\n";
          ++parameter_error_count_;
          value = TIME_CONSTANT_MINIMUM;
          return false;
        }
        if (value < ZERO<RealT>)
        {
          Log::error() << "Reecb: " << name << " must be non-negative\n";
          ++parameter_error_count_;
          value = TIME_CONSTANT_MINIMUM;
          return false;
        }

        const bool raised = value < TIME_CONSTANT_MINIMUM;
        value             = std::max(value, TIME_CONSTANT_MINIMUM);
        return raised;
      }

      /**
       * @brief Read parameters from model data
       *
       * Omitted parameters retain their documented defaults. Loading errors
       * are counted for verify() rather than thrown.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;
        mva_given_             = data.parameters.contains(Params::mva);
        Vref0_given_           = false;

        loadRealParameter(data, Params::mva, mva_base_, "mva");
        loadBooleanParameter(data, Params::PfFlag, PfFlag_, "PfFlag");
        loadBooleanParameter(data, Params::VFlag, VFlag_, "VFlag");
        loadBooleanParameter(data, Params::QFlag, QFlag_, "QFlag");
        loadBooleanParameter(data, Params::Pqflag, Pqflag_, "Pqflag");
        loadRealParameter(data, Params::Trv, Trv_, "Trv");
        loadRealParameter(data, Params::Tp, Tp_, "Tp");
        if (data.parameters.contains(Params::Vref0))
        {
          loadRealParameter(data, Params::Vref0, Vref0_, "Vref0");
          Vref0_given_ = true;
        }
        loadRealParameter(data, Params::Vdip, Vdip_, "Vdip");
        loadRealParameter(data, Params::Vup, Vup_, "Vup");
        loadRealParameter(data, Params::dbd1, dbd1_, "dbd1");
        loadRealParameter(data, Params::dbd2, dbd2_, "dbd2");
        loadRealParameter(data, Params::kqv, kqv_, "kqv");
        loadRealParameter(data, Params::Iql1, Iql1_, "Iql1");
        loadRealParameter(data, Params::Iqh1, Iqh1_, "Iqh1");
        loadRealParameter(data, Params::Qmax, Qmax_, "Qmax");
        loadRealParameter(data, Params::Qmin, Qmin_, "Qmin");
        loadRealParameter(data, Params::Kqp, Kqp_, "Kqp");
        loadRealParameter(data, Params::Kqi, Kqi_, "Kqi");
        loadRealParameter(data, Params::Vmax, Vmax_, "Vmax");
        loadRealParameter(data, Params::Vmin, Vmin_, "Vmin");
        loadRealParameter(data, Params::Kvp, Kvp_, "Kvp");
        loadRealParameter(data, Params::Kvi, Kvi_, "Kvi");
        loadRealParameter(data, Params::Tiq, Tiq_, "Tiq");
        loadRealParameter(data, Params::Tpord, Tpord_, "Tpord");
        loadRealParameter(data, Params::dPmax, dPmax_, "dPmax");
        loadRealParameter(data, Params::dPmin, dPmin_, "dPmin");
        loadRealParameter(data, Params::Pmax, Pmax_, "Pmax");
        loadRealParameter(data, Params::Pmin, Pmin_, "Pmin");
        loadRealParameter(data, Params::Imax, Imax_, "Imax");

        setDerivedParameters();
      }

      /**
       * @brief Bind monitor selections to REECB internal states
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::iqcmd, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::IQCMD)]; });
        monitor_->set(Variable::ipcmd, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::IPCMD)]; });
        monitor_->set(Variable::vmeas, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::VMEAS)]; });
        monitor_->set(Variable::pmeas, [this]
                      { return y_.getData()[static_cast<size_t>(ReecbInternalVariables::PMEAS)]; });
      }

      /**
       * @brief Resolve parameter-derived constants and selector masks
       *
       * Raises explicit controller lags in place, converts any supplied component
       * rating, and resolves selector masks. Invalid lag inputs are
       * recorded before replacement so verify() retains each error.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::setDerivedParameters()
      {
        bool floor_warning = false;

        floor_warning |= floorTimeConstant(Trv_, "Trv");
        floor_warning |= floorTimeConstant(Tp_, "Tp");
        floor_warning |= floorTimeConstant(Tiq_, "Tiq");
        floor_warning |= floorTimeConstant(Tpord_, "Tpord");

        if (floor_warning)
        {
          Log::warning() << "Reecb: any of Trv, Tp, Tiq, or Tpord below "
                         << TIME_CONSTANT_MINIMUM
                         << " s is raised to that floor to keep the controller lags well posed\n";
        }

        va_component_base_ = mva_base_ * static_cast<RealT>(1.0e6);

        pf_on_ = ZERO<RealT>;
        if (PfFlag_)
        {
          pf_on_ = ONE<RealT>;
        }
        pf_off_ = ONE<RealT> - pf_on_;

        q_on_ = ZERO<RealT>;
        if (QFlag_)
        {
          q_on_ = ONE<RealT>;
        }
        q_off_ = ONE<RealT> - q_on_;

        q_pi_on_ = ZERO<RealT>;
        if (QFlag_ && VFlag_)
        {
          q_pi_on_ = ONE<RealT>;
        }

        v_ref_on_ = ZERO<RealT>;
        if (QFlag_ && !VFlag_)
        {
          v_ref_on_ = ONE<RealT>;
        }
        q_ref_on_ = ONE<RealT> - v_ref_on_;

        pq_on_ = ZERO<RealT>;
        if (Pqflag_)
        {
          pq_on_ = ONE<RealT>;
        }
        pq_off_ = ONE<RealT> - pq_on_;
      }

      /**
       * @brief Evaluate log(1 - exp(-x)) accurately for positive x
       *
       * The small-x hyperbolic form avoids cancellation; the large-x form uses
       * log1p. The algebraically equivalent branches agree in value and first
       * derivative at x = log(2).
       *
       * @param[in] x Strictly positive argument.
       * @return Numerically stable value of log(1 - exp(-x)).
       * @pre `x > 0`.
       * @warning This function contains conditional branching and as such can
       * be used in initialization methods but not in residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::logOneMinusExp(RealT x) const
      {
        static constexpr RealT log_two = std::numbers::ln2_v<RealT>;

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      /**
       * @brief Recover the input that produces a requested smooth-clamp output
       *
       * @param[in] output Requested output.
       * @param[in] lower Lower smooth-clamp limit.
       * @param[in] upper Upper smooth-clamp limit.
       * @return Inverse of `Math::clamp`, or a nonfinite value when no input
       *         produces `output`.
       * @pre `lower <= upper`.
       * @warning This function contains conditional branching and as such can
       * be used in initialization methods but not in residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::unclamp(RealT output, RealT lower, RealT upper) const
      {
        const RealT a = Math::MU<RealT> * (output - lower);
        const RealT b = Math::MU<RealT> * (upper - output);
        return lower + (a + logOneMinusExp(a) - logOneMinusExp(b)) / Math::MU<RealT>;
      }

      /**
       * @brief Resolve the REECB component power base
       *
       * @return The supplied component base, or the system base when `mva` is omitted.
       */
      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::componentPowerBase() const
      {
        if (mva_given_)
        {
          return va_component_base_;
        }
        return va_system_base_;
      }

      /**
       * @brief Convert a system-base power or current to REECB component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the REECB component base.
       */
      template <typename scalar_type, typename index_type>
      template <typename ValueT>
      [[gnu::always_inline]] inline ValueT
      Reecb<scalar_type, index_type>::toComponentBase(ValueT value) const
      {
        return value * (va_system_base_ / componentPowerBase());
      }

      /**
       * @brief Convert a component-base power or current to system base
       *
       * @param[in] value Quantity on the REECB component base.
       * @return The same quantity on the system base.
       */
      template <typename scalar_type, typename index_type>
      template <typename ValueT>
      ValueT Reecb<scalar_type, index_type>::toSystemBase(ValueT value) const
      {
        return value * (componentPowerBase() / va_system_base_);
      }

      /**
       * @brief Access the terminal-bus real voltage component
       *
       * @return Mutable reference to the bus real voltage state.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Reecb<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /**
       * @brief Access the terminal-bus imaginary voltage component
       *
       * @return Mutable reference to the bus imaginary voltage state.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Reecb<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
