/**
 * @file ReecbImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REECB electrical-control model.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <mutex>
#include <numbers>
#include <numeric>
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

      template <typename scalar_type, typename index_type>
      struct Reecb<scalar_type, index_type>::InitialPoint
      {
        std::array<RealT, static_cast<size_t>(ReecbInternalVariables::IQCMD)>   variables{};
        std::array<RealT, static_cast<size_t>(ReecbExternalVariables::MAXIMUM)> signals{};
        RealT                                                                   qmin{};
        RealT                                                                   qmax{};
        RealT                                                                   vmin{};
        RealT                                                                   vmax{};
        RealT                                                                   pmin{};
        RealT                                                                   pmax{};
        RealT                                                                   imax{};
        RealT                                                                   vref{};
      };

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

        this->allocateExternalVectors(static_cast<IdxT>(ReecbExternalVariables::MAXIMUM));

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

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Reecb: " << message << '\n';
            ret += 1;
          }
        };

        check(bus_ != nullptr, "terminal bus is required");

        const RealT component_power_base = componentPowerBase();
        const bool  valid_component_base = std::isfinite(component_power_base) && component_power_base > ZERO<RealT>;
        const bool  valid_system_base    = std::isfinite(va_system_base_) && va_system_base_ > ZERO<RealT>;
        check(valid_component_base, "component power base must be finite and positive");
        check(valid_system_base, "system power base must be finite and positive");
        if (valid_component_base && valid_system_base)
        {
          const RealT system_to_component = va_system_base_ / component_power_base;
          const RealT component_to_system = component_power_base / va_system_base_;
          check(
              std::isfinite(system_to_component)
                  && system_to_component > ZERO<RealT>
                  && std::isfinite(component_to_system)
                  && component_to_system > ZERO<RealT>,
              "system/component power-base conversion ratios must be finite and positive");
        }

        check(std::isfinite(Trv_), "Trv must be finite");
        check(std::isfinite(Tp_), "Tp must be finite");
        check(std::isfinite(Vref0_), "Vref0 must be finite");

        const bool finite_voltage_thresholds = std::isfinite(Vdip_) && std::isfinite(Vup_);
        check(finite_voltage_thresholds, "Vdip and Vup must be finite");
        if (finite_voltage_thresholds)
        {
          check(Vdip_ < Vup_, "Vdip must be less than Vup");
        }

        const bool finite_voltage_deadband = std::isfinite(dbd1_) && std::isfinite(dbd2_);
        check(finite_voltage_deadband, "dbd1 and dbd2 must be finite");
        if (finite_voltage_deadband)
        {
          check(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required");
        }

        check(std::isfinite(kqv_) && kqv_ >= ZERO<RealT>, "kqv must be finite and non-negative");

        const bool finite_injection_limits = std::isfinite(Iql1_) && std::isfinite(Iqh1_);
        check(finite_injection_limits, "Iql1 and Iqh1 must be finite");
        if (finite_injection_limits)
        {
          check(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1");
        }

        const bool finite_reactive_limits = std::isfinite(Qmin_) && std::isfinite(Qmax_);
        check(finite_reactive_limits, "Qmin and Qmax must be finite");
        if (finite_reactive_limits)
        {
          check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        }

        check(std::isfinite(Kqp_) && Kqp_ >= ZERO<RealT>, "Kqp must be finite and non-negative");
        check(std::isfinite(Kqi_) && Kqi_ >= ZERO<RealT>, "Kqi must be finite and non-negative");

        const bool finite_voltage_limits = std::isfinite(Vmin_) && std::isfinite(Vmax_);
        check(finite_voltage_limits, "Vmin and Vmax must be finite");
        if (finite_voltage_limits)
        {
          check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        }

        check(std::isfinite(Kvp_) && Kvp_ >= ZERO<RealT>, "Kvp must be finite and non-negative");
        check(std::isfinite(Kvi_) && Kvi_ >= ZERO<RealT>, "Kvi must be finite and non-negative");
        check(std::isfinite(Tiq_), "Tiq must be finite");
        check(std::isfinite(Tpord_), "Tpord must be finite");

        const bool finite_ramp_limits = std::isfinite(dPmin_) && std::isfinite(dPmax_);
        check(finite_ramp_limits, "dPmin and dPmax must be finite");
        if (finite_ramp_limits)
        {
          check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
        }

        const bool finite_active_limits = std::isfinite(Pmin_) && std::isfinite(Pmax_);
        check(finite_active_limits, "Pmin and Pmax must be finite");
        if (finite_active_limits)
        {
          check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
        }

        check(std::isfinite(Imax_) && Imax_ > ZERO<RealT>, "Imax must be finite and positive");

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
       * power feedback or reconstructs unattached feedback, and constructs the
       * remaining states and reference setpoints at a steady operating point.
       *
       * @pre allocate() has completed.
       * @pre verify() reports a valid parameter and port configuration.
       * @pre The terminal bus and current-command states are initialized.
       *
       * @post On failure no state, derivative, parameter, or signal storage is
       *       modified.
       *
       * @return 0 on success; nonzero when allocation, configuration,
       *         initial-value, limiter inversion, or steady-state checks fail.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initialize()
      {
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

        InitialPoint point;
        if (!buildInitialPoint(point))
        {
          return 1;
        }

        commitInitialPoint(point);
        return 0;
      }

      /**
       * @brief Construct and validate a steady initial point
       *
       * @param[out] point Initial point.
       * @return true on success.
       */
      template <typename scalar_type, typename index_type>
      bool Reecb<scalar_type, index_type>::buildInitialPoint(InitialPoint& point)
      {
        const auto VMEAS  = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS  = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ   = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV   = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV     = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD   = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT     = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto VSAFE  = static_cast<size_t>(ReecbInternalVariables::VSAFE);
        const auto SDIP   = static_cast<size_t>(ReecbInternalVariables::SDIP);
        const auto IQV    = static_cast<size_t>(ReecbInternalVariables::IQV);
        const auto QREF   = static_cast<size_t>(ReecbInternalVariables::QREF);
        const auto EQ     = static_cast<size_t>(ReecbInternalVariables::EQ);
        const auto VPIQ   = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        const auto EPIV   = static_cast<size_t>(ReecbInternalVariables::EPIV);
        const auto RPORD  = static_cast<size_t>(ReecbInternalVariables::RPORD);
        const auto ILMAX  = static_cast<size_t>(ReecbInternalVariables::ILMAX);
        const auto ILCAP  = static_cast<size_t>(ReecbInternalVariables::ILCAP);
        const auto IQMAX  = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        const auto IPMAX  = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        const auto IQBASE = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        const auto IQRAW  = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        const auto IQCMD  = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD  = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        const auto* y = y_.getData();

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
          return false;
        }
        if (vt0 <= ZERO<RealT>)
        {
          Log::error() << "Reecb: initial terminal-voltage magnitude must be positive\n";
          return false;
        }

        if (ipcmd0 < ZERO<RealT>)
        {
          Log::error() << "Reecb: initial active-current command must be non-negative\n";
          return false;
        }

        const RealT verr0   = Math::deadband2(vref0 - vmeas0, dbd1_, dbd2_);
        const RealT iqv0    = Math::clamp(kqv_ * verr0, Iql1_, Iqh1_);
        const RealT iqabs0  = std::abs(iqcmd0);
        RealT       iqneed0 = iqabs0;
        if (QFlag_ && iqabs0 > ZERO<RealT>)
        {
          iqneed0 += std::numbers::ln2_v<RealT> / Math::MU<RealT> + INITIALIZATION_TOLERANCE;
        }

        RealT high0 = iqabs0;
        RealT low0  = ipcmd0;
        if (Pqflag_)
        {
          high0 = ipcmd0;
          low0  = iqneed0;
        }
        // Q priority uses Imax directly for reactive current, so include the
        // smooth-clamp recovery margin carried by iqneed0.
        const auto current_limit = solveInitialLimit(
            std::max({Imax_, high0, low0, iqneed0}), high0, low0);
        if (!current_limit)
        {
          Log::error() << "Reecb: adjusted Imax cannot include the initial current commands\n";
          return false;
        }
        const RealT imax   = current_limit->total_limit;
        const RealT ilmax0 = current_limit->continuation;
        const RealT ilcap0 = current_limit->off_axis_capacity;

        RealT iqmax0 = imax;
        RealT ipmax0 = ilcap0;
        if (Pqflag_)
        {
          iqmax0 = ilcap0;
          ipmax0 = imax;
        }

        RealT ipraw0 = ZERO<RealT>;
        RealT iqraw0 = ZERO<RealT>;
        if (!iclamp(ipcmd0, ZERO<RealT>, ipmax0, ipraw0)
            || !iclamp(iqcmd0, -iqmax0, iqmax0, iqraw0))
        {
          Log::error() << "Reecb: initial current commands cannot be reproduced by their limiters\n";
          return false;
        }

        const RealT iqctl0 = iqraw0 - iqv0;
        const RealT pord0  = vmeas_safe0 * ipraw0;
        RealT       qmin   = Qmin_;
        RealT       qmax   = Qmax_;
        RealT       vmin   = Vmin_;
        RealT       vmax   = Vmax_;
        if (QFlag_ && VFlag_)
        {
          qmin = std::min(Qmin_, qgen0);
          qmax = std::max(Qmax_, qgen0);
          vmin = std::min(Vmin_, vmeas0);
          vmax = std::max(Vmax_, vmeas0);

          const RealT infinity = std::numeric_limits<RealT>::infinity();
          if (qmin == qgen0 && qmin < qmax)
          {
            qmin = std::nextafter(qmin, -infinity);
          }
          if (qmax == qgen0 && qmin < qmax)
          {
            qmax = std::nextafter(qmax, infinity);
          }
          if (vmin == vmeas0 && vmin < vmax)
          {
            vmin = std::nextafter(vmin, -infinity);
          }
          if (vmax == vmeas0 && vmin < vmax)
          {
            vmax = std::nextafter(vmax, infinity);
          }
        }
        const RealT pmin         = std::min(Pmin_, pord0);
        const RealT pmax         = std::max(Pmax_, pord0);
        const RealT pref0_system = toSystemBase(pord0);

        RealT qtarget0 = ZERO<RealT>;
        if (!QFlag_)
        {
          qtarget0 = iqctl0 * vmeas_safe0;
        }
        else if (VFlag_ && !iclamp(qgen0, qmin, qmax, qtarget0))
        {
          Log::error() << "Reecb: reactive-power limiter has no finite steady input\n";
          return false;
        }

        RealT qref0      = ZERO<RealT>;
        RealT qext0_port = ZERO<RealT>;
        RealT pfaref0    = ZERO<RealT>;

        if (QFlag_ && !VFlag_)
        {
          qext0_port = vmeas0;
        }
        else if (PfFlag_)
        {
          if (pmeas0 == ZERO<RealT> && qtarget0 != ZERO<RealT>)
          {
            Log::error() << "Reecb: power-factor mode cannot reproduce the reactive target at zero active power\n";
            return false;
          }
          if (pmeas0 != ZERO<RealT>)
          {
            pfaref0 = std::atan(qtarget0 / pmeas0);
          }
          qref0 = pmeas0 * std::tan(pfaref0);
          if (std::abs(qref0 - qtarget0) > std::abs(qtarget0) * INITIALIZATION_TOLERANCE)
          {
            Log::error() << "Reecb: power-factor angle cannot reproduce the reactive target\n";
            return false;
          }
        }
        else
        {
          qext0_port = toSystemBase(qtarget0);
          qref0      = toComponentBase(qext0_port);
        }

        const RealT eq0   = Math::clamp(qref0, qmin, qmax) - qgen0;
        RealT       xpiq0 = ZERO<RealT>;
        if (QFlag_ && VFlag_)
        {
          RealT vpiq_input0 = ZERO<RealT>;
          if (!iclamp(vmeas0, vmin, vmax, vpiq_input0))
          {
            Log::error() << "Reecb: voltage limiter has no finite steady input\n";
            return false;
          }
          xpiq0 = vpiq_input0 - Kqp_ * eq0;
        }

        const RealT vpiq0 = Math::clamp(Kqp_ * eq0 + xpiq0, vmin, vmax);
        RealT       epiv0 = ZERO<RealT>;
        RealT       qv0   = ZERO<RealT>;
        RealT       xpiv0 = ZERO<RealT>;

        if (QFlag_)
        {
          if (VFlag_)
          {
            epiv0 = vpiq0 - vmeas0;
          }
          else
          {
            epiv0 = qext0_port - vmeas0;
          }

          if (iqmax0 <= INITIALIZATION_TOLERANCE)
          {
            xpiv0 = -Kvp_ * epiv0;
          }
          else
          {
            RealT iqctl_input0 = ZERO<RealT>;
            if (!iclamp(iqctl0, -iqmax0, iqmax0, iqctl_input0))
            {
              Log::error() << "Reecb: voltage-controller current cannot be reproduced by its limiter\n";
              return false;
            }
            xpiv0 = iqctl_input0 - Kvp_ * epiv0;
          }
        }
        else
        {
          qv0 = qref0 / vmeas_safe0;
        }

        const RealT sdip0  = Math::inside(vt0, Vdip_, Vup_);
        RealT       qrate0 = ZERO<RealT>;
        if (QFlag_ && VFlag_)
        {
          qrate0 = sdip0 * Math::antiwindup(Kqp_ * eq0 + xpiq0, Kqi_ * eq0, vmin, vmax);
        }

        const ScalarT vstate0{Kvp_ * epiv0 + xpiv0};
        const ScalarT vderiv0{Kvi_ * epiv0};
        RealT         vrate0 = ZERO<RealT>;
        if (QFlag_)
        {
          vrate0 = sdip0 * static_cast<RealT>(awband(vstate0, vderiv0, ScalarT{iqmax0}));
        }

        const RealT iqbase0     = Math::clamp(Kvp_ * epiv0 + xpiv0, -iqmax0, iqmax0);
        RealT       iqraw_check = qv0 + iqv0;
        if (QFlag_)
        {
          iqraw_check = iqbase0 + iqv0;
        }

        const RealT iqcmd_check = Math::clamp(iqraw_check, -iqmax0, iqmax0);
        const RealT ipcmd_check = Math::clamp(pord0 / vmeas_safe0, ZERO<RealT>, ipmax0);

        if (!std::isfinite(imax) || !std::isfinite(ilmax0) || !std::isfinite(ilcap0)
            || !std::isfinite(iqmax0) || !std::isfinite(ipmax0)
            || !std::isfinite(ipraw0) || !std::isfinite(iqraw0) || !std::isfinite(pord0)
            || !std::isfinite(pref0_system) || !std::isfinite(qtarget0) || !std::isfinite(qref0)
            || !std::isfinite(qext0_port) || !std::isfinite(pfaref0) || !std::isfinite(eq0)
            || !std::isfinite(xpiq0) || !std::isfinite(epiv0) || !std::isfinite(xpiv0)
            || !std::isfinite(qv0) || !std::isfinite(qrate0) || !std::isfinite(vrate0)
            || !std::isfinite(iqbase0) || !std::isfinite(iqraw_check)
            || !std::isfinite(iqcmd_check) || !std::isfinite(ipcmd_check))
        {
          Log::error() << "Reecb: initialization produced a nonfinite value\n";
          return false;
        }
        if (std::abs(qrate0) > INITIALIZATION_TOLERANCE || std::abs(vrate0) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Reecb: controller state rate is nonzero at initialization\n";
          return false;
        }
        if (std::abs(iqcmd_check - iqcmd0) > INITIALIZATION_TOLERANCE
            || std::abs(ipcmd_check - ipcmd0) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Reecb: current-command limiter reconstruction is inexact\n";
          return false;
        }

        point.variables[VMEAS]  = vmeas0;
        point.variables[PMEAS]  = pmeas0;
        point.variables[XPIQ]   = xpiq0;
        point.variables[XPIV]   = xpiv0;
        point.variables[QV]     = qv0;
        point.variables[PORD]   = pord0;
        point.variables[VT]     = vt0;
        point.variables[VSAFE]  = vmeas_safe0;
        point.variables[SDIP]   = sdip0;
        point.variables[IQV]    = iqv0;
        point.variables[QREF]   = qref0;
        point.variables[EQ]     = eq0;
        point.variables[VPIQ]   = vpiq0;
        point.variables[EPIV]   = epiv0;
        point.variables[RPORD]  = ZERO<RealT>;
        point.variables[ILMAX]  = ilmax0;
        point.variables[ILCAP]  = ilcap0;
        point.variables[IQMAX]  = iqmax0;
        point.variables[IPMAX]  = ipmax0;
        point.variables[IQBASE] = iqbase0;
        point.variables[IQRAW]  = iqraw_check;
        point.signals[PE]       = pe0_system;
        point.signals[QGEN]     = qgen0_system;
        point.signals[QEXT]     = qext0_port;
        point.signals[PFAREF]   = pfaref0;
        point.signals[PREF]     = pref0_system;

        point.qmin = qmin;
        point.qmax = qmax;
        point.vmin = vmin;
        point.vmax = vmax;
        point.pmin = pmin;
        point.pmax = pmax;
        point.imax = imax;
        point.vref = vref0;
        return true;
      }

      /**
       * @brief Commit a validated initial point
       *
       * @param[in] point Initial point.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::commitInitialPoint(const InitialPoint& point)
      {
        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        const bool q_adjusted    = point.qmin != Qmin_ || point.qmax != Qmax_;
        const bool v_adjusted    = point.vmin != Vmin_ || point.vmax != Vmax_;
        const bool p_adjusted    = point.pmin != Pmin_ || point.pmax != Pmax_;
        const bool imax_adjusted = point.imax != Imax_;

        Qmin_ = point.qmin;
        Qmax_ = point.qmax;
        Vmin_ = point.vmin;
        Vmax_ = point.vmax;
        Pmin_ = point.pmin;
        Pmax_ = point.pmax;
        Imax_ = point.imax;

        auto* y = y_.getData();
        // Preserve the owned IQCMD and IPCMD inputs.
        for (size_t entry = 0; entry < point.variables.size(); ++entry)
        {
          y[entry] = static_cast<ScalarT>(point.variables[entry]);
        }

        if (!Vref0_given_)
        {
          Vref0_ = point.vref;
        }

        pe_set_     = static_cast<ScalarT>(point.signals[PE]);
        qgen_set_   = static_cast<ScalarT>(point.signals[QGEN]);
        qext_set_   = static_cast<ScalarT>(point.signals[QEXT]);
        pfaref_set_ = static_cast<ScalarT>(point.signals[PFAREF]);
        pref_set_   = static_cast<ScalarT>(point.signals[PREF]);

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

        if (q_adjusted)
        {
          Log::warning() << "Reecb: Qmin/Qmax adjusted to include the initial reactive power\n";
        }
        if (v_adjusted)
        {
          Log::warning() << "Reecb: Vmin/Vmax adjusted to include the initial terminal voltage\n";
        }
        if (p_adjusted)
        {
          Log::warning() << "Reecb: Pmin/Pmax adjusted to include the initial active-power order\n";
        }
        if (imax_adjusted)
        {
          Log::warning() << "Reecb: Imax adjusted to include the initial current commands\n";
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
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
       * @brief Gather external variables and index maps.
       *
       * Unattached signal inputs retain the values latched by initialize().
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::gatherExternalVariables()
      {
        auto* y_ext = y_ext_.getData();

        const auto VR_EXT = static_cast<size_t>(ReecbExternalVariables::VR);
        const auto VI_EXT = static_cast<size_t>(ReecbExternalVariables::VI);
        const auto PE     = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN   = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto QEXT   = static_cast<size_t>(ReecbExternalVariables::QEXT);
        const auto PFAREF = static_cast<size_t>(ReecbExternalVariables::PFAREF);
        const auto PREF   = static_cast<size_t>(ReecbExternalVariables::PREF);

        y_ext[VR_EXT]                = Vr();
        y_ext[VI_EXT]                = Vi();
        variable_indices_ext_[VR_EXT] = INVALID_INDEX<IdxT>;
        variable_indices_ext_[VI_EXT] = INVALID_INDEX<IdxT>;
        if (bus_->size() > 0)
        {
          variable_indices_ext_[VR_EXT] = bus_->getVariableIndex(0);
          variable_indices_ext_[VI_EXT] = bus_->getVariableIndex(1);
        }

        y_ext[PE]                    = pe_set_;
        y_ext[QGEN]                  = qgen_set_;
        y_ext[QEXT]                  = qext_set_;
        y_ext[PFAREF]                = pfaref_set_;
        y_ext[PREF]                  = pref_set_;
        variable_indices_ext_[PE]     = INVALID_INDEX<IdxT>;
        variable_indices_ext_[QGEN]   = INVALID_INDEX<IdxT>;
        variable_indices_ext_[QEXT]   = INVALID_INDEX<IdxT>;
        variable_indices_ext_[PFAREF] = INVALID_INDEX<IdxT>;
        variable_indices_ext_[PREF]   = INVALID_INDEX<IdxT>;

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          y_ext[PE] = signals_.template readExternalVariable<ReecbExternalVariables::PE>();
          variable_indices_ext_[PE] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PE>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          y_ext[QGEN] = signals_.template readExternalVariable<ReecbExternalVariables::QGEN>();
          variable_indices_ext_[QGEN] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::QGEN>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          y_ext[QEXT] = signals_.template readExternalVariable<ReecbExternalVariables::QEXT>();
          variable_indices_ext_[QEXT] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::QEXT>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          y_ext[PFAREF] =
              signals_.template readExternalVariable<ReecbExternalVariables::PFAREF>();
          variable_indices_ext_[PFAREF] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PFAREF>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          y_ext[PREF] = signals_.template readExternalVariable<ReecbExternalVariables::PREF>();
          variable_indices_ext_[PREF] =
              signals_.template readExternalVariableIndex<ReecbExternalVariables::PREF>();
        }
      }

      /**
       * @brief Evaluate the internal REECB residual equations.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateInternalResidual()
      {
        gatherExternalVariables();

        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.getData(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Evaluate internal equations and external contributions.
       *
       * REECB contributes no external residual, so the base implementation
       * returns zero after the internal equations are evaluated.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateResidual()
      {
        evaluateInternalResidual();
        return this->evaluateExternalResidual();
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
       * Evaluates the residual rows in enum order. Selector masks keep the
       * differentiated path branch-free.
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] y_ext External variable values.
       * @param[out] f Internal residuals.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline int
      Reecb<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* y_ext,
          ScalarT*       f)
      {
        const auto VMEAS  = static_cast<size_t>(ReecbInternalVariables::VMEAS);
        const auto PMEAS  = static_cast<size_t>(ReecbInternalVariables::PMEAS);
        const auto XPIQ   = static_cast<size_t>(ReecbInternalVariables::XPIQ);
        const auto XPIV   = static_cast<size_t>(ReecbInternalVariables::XPIV);
        const auto QV     = static_cast<size_t>(ReecbInternalVariables::QV);
        const auto PORD   = static_cast<size_t>(ReecbInternalVariables::PORD);
        const auto VT     = static_cast<size_t>(ReecbInternalVariables::VT);
        const auto VSAFE  = static_cast<size_t>(ReecbInternalVariables::VSAFE);
        const auto SDIP   = static_cast<size_t>(ReecbInternalVariables::SDIP);
        const auto IQV    = static_cast<size_t>(ReecbInternalVariables::IQV);
        const auto QREF   = static_cast<size_t>(ReecbInternalVariables::QREF);
        const auto EQ     = static_cast<size_t>(ReecbInternalVariables::EQ);
        const auto VPIQ   = static_cast<size_t>(ReecbInternalVariables::VPIQ);
        const auto EPIV   = static_cast<size_t>(ReecbInternalVariables::EPIV);
        const auto RPORD  = static_cast<size_t>(ReecbInternalVariables::RPORD);
        const auto ILMAX  = static_cast<size_t>(ReecbInternalVariables::ILMAX);
        const auto ILCAP  = static_cast<size_t>(ReecbInternalVariables::ILCAP);
        const auto IQMAX  = static_cast<size_t>(ReecbInternalVariables::IQMAX);
        const auto IPMAX  = static_cast<size_t>(ReecbInternalVariables::IPMAX);
        const auto IQBASE = static_cast<size_t>(ReecbInternalVariables::IQBASE);
        const auto IQRAW  = static_cast<size_t>(ReecbInternalVariables::IQRAW);
        const auto IQCMD  = static_cast<size_t>(ReecbInternalVariables::IQCMD);
        const auto IPCMD  = static_cast<size_t>(ReecbInternalVariables::IPCMD);

        const auto VR_EXT = static_cast<size_t>(ReecbExternalVariables::VR);
        const auto VI_EXT = static_cast<size_t>(ReecbExternalVariables::VI);
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
        const ScalarT vsafe        = y[VSAFE];
        const ScalarT sdip         = y[SDIP];
        const ScalarT iqv          = y[IQV];
        const ScalarT qref         = y[QREF];
        const ScalarT eq           = y[EQ];
        const ScalarT vpiq         = y[VPIQ];
        const ScalarT epiv         = y[EPIV];
        const ScalarT rpord        = y[RPORD];
        const ScalarT ilmax        = y[ILMAX];
        const ScalarT ilcap        = y[ILCAP];
        const ScalarT iqmax        = y[IQMAX];
        const ScalarT ipmax        = y[IPMAX];
        const ScalarT iqbase       = y[IQBASE];
        const ScalarT iqraw        = y[IQRAW];
        const ScalarT iqcmd_system = y[IQCMD];
        const ScalarT ipcmd_system = y[IPCMD];

        const ScalarT vmeas_dot = yp[VMEAS];
        const ScalarT pmeas_dot = yp[PMEAS];
        const ScalarT xpiq_dot  = yp[XPIQ];
        const ScalarT xpiv_dot  = yp[XPIV];
        const ScalarT qv_dot    = yp[QV];
        const ScalarT pord_dot  = yp[PORD];

        const ScalarT vr = y_ext[VR_EXT];
        const ScalarT vi = y_ext[VI_EXT];

        const ScalarT pe     = toComponentBase(y_ext[PE]);
        const ScalarT qgen   = toComponentBase(y_ext[QGEN]);
        const ScalarT extref = y_ext[QEXT];
        const ScalarT pfaref = y_ext[PFAREF];
        const ScalarT pref   = toComponentBase(y_ext[PREF]);
        const ScalarT iqcmd  = toComponentBase(iqcmd_system);
        const ScalarT ipcmd  = toComponentBase(ipcmd_system);

        const ScalarT verr        = Math::deadband2(Vref0_ - vmeas, dbd1_, dbd2_);
        const ScalarT q_pi_state  = Kqp_ * eq + xpiq;
        const ScalarT v_pi_state  = Kvp_ * epiv + xpiv;
        const ScalarT fpord       = (pref - pord) / Tpord_;
        const ScalarT ilnorm      = std::sqrt(ilmax * ilmax + INITIALIZATION_TOLERANCE);
        // Select before the factored square to avoid 0 * inf on the inactive path.
        const ScalarT high        = pq_on_ * ipcmd + pq_off_ * iqcmd;
        const ScalarT q_pi_rate   = q_pi_on_ * sdip * Math::antiwindup(q_pi_state, Kqi_ * eq, Vmin_, Vmax_);
        const ScalarT v_pi_rate   = q_on_ * sdip * awband(v_pi_state, Kvi_ * epiv, iqmax);
        const ScalarT qv_rate     = q_off_ * sdip * (qref / vsafe - qv) / Tiq_;
        const ScalarT pord_rate   = sdip * Math::antiwindup(pord, rpord, Pmin_, Pmax_);
        const ScalarT iqv_target  = Math::clamp(kqv_ * verr, Iql1_, Iqh1_);
        // The Volt/VAr channel is a system-base reactive power unless
        // direct-voltage mode selects it as a terminal-voltage reference,
        // which takes no power-base conversion.
        const ScalarT qref_target = q_ref_on_ * (pf_on_ * pmeas * std::tan(pfaref) + pf_off_ * toComponentBase(extref));

        f[VMEAS]  = -vmeas_dot + (vt - vmeas) / Trv_;
        f[PMEAS]  = -pmeas_dot + (pe - pmeas) / Tp_;
        f[XPIQ]   = -xpiq_dot + q_pi_rate;
        f[XPIV]   = -xpiv_dot + v_pi_rate;
        f[QV]     = -qv_dot + qv_rate;
        f[PORD]   = -pord_dot + pord_rate;
        f[VT]     = -vt * vt + vr * vr + vi * vi;
        f[VSAFE]  = -vsafe + Math::max(vmeas, VMEAS_MINIMUM);
        f[SDIP]   = -sdip + Math::inside(vt, Vdip_, Vup_);
        f[IQV]    = -iqv + iqv_target;
        f[QREF]   = -qref + qref_target;
        f[EQ]     = -eq + Math::clamp(qref, Qmin_, Qmax_) - qgen;
        f[VPIQ]   = -vpiq + Math::clamp(q_pi_state, Vmin_, Vmax_);
        f[EPIV]   = -epiv + q_pi_on_ * vpiq + v_ref_on_ * extref - q_on_ * vmeas;
        f[RPORD]  = -rpord + aslew(fpord, dPmin_, dPmax_);
        f[ILMAX]  = -ilmax * ilnorm + (Imax_ - high) * (Imax_ + high);
        f[ILCAP]  = -ilcap + (ilmax / ilnorm) * ilmax;
        f[IQMAX]  = -iqmax + pq_on_ * ilcap + pq_off_ * Imax_;
        f[IPMAX]  = -ipmax + pq_on_ * Imax_ + pq_off_ * ilcap;
        f[IQBASE] = -iqbase + Math::clamp(v_pi_state, -iqmax, iqmax);
        f[IQRAW]  = -iqraw + q_on_ * iqbase + q_off_ * qv + iqv;
        f[IQCMD]  = -iqcmd + Math::clamp(iqraw, -iqmax, iqmax);
        f[IPCMD]  = -ipcmd + Math::clamp(pord / vsafe, ZERO<RealT>, ipmax);

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Smooth asymmetric slew-rate limiter
       *
       * @param[in] rate Unconstrained rate.
       * @param[in] lower Negative rate limit.
       * @param[in] upper Positive rate limit.
       * @return Limited rate.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Reecb<scalar_type, index_type>::aslew(ScalarT rate, RealT lower, RealT upper)
      {
        assert(lower < ZERO<RealT> && ZERO<RealT> < upper);
        return rate
               / (ONE<RealT> + Math::ramp(rate / upper - ONE<RealT>) + Math::ramp(rate / lower - ONE<RealT>));
      }

      /**
       * @brief Smooth anti-windup derivative within a moving symmetric band
       *
       * Math::antiwindup over [-band, band] with a band edge that is an
       * algebraic quantity, so differentiation carries the band's own
       * contributions through the gate.
       *
       * @param[in] state Limited PI state.
       * @param[in] rate Pre-limit derivative of state.
       * @param[in] band Nonnegative symmetric band edge.
       * @return Anti-windup-limited derivative.
       *
       * @todo Fold moving-limit support into Math::antiwindup in CommonMath.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      Reecb<scalar_type, index_type>::awband(ScalarT state, ScalarT rate, ScalarT band)
      {
        const ScalarT above_min = Math::above(state, -band);
        const ScalarT below_max = Math::below(state, band);
        return (above_min * below_max + (ONE<RealT> - below_max) * Math::sigmoid(-rate)
                + (ONE<RealT> - above_min) * Math::sigmoid(rate))
               * rate;
      }

      /**
       * @brief Compute the initial current-circle continuation state
       *
       * Solves the implemented `ILMAX` row for its nonnegative continuation
       * state at a total-current limit and priority-axis current.
       *
       * @param[in] imax Total-current limit on the component base.
       * @param[in] high Priority-axis current command on the component base.
       * @return Initial `ILMAX` state on the component base, or a quiet NaN
       *         when the circle geometry is invalid or unrepresentable.
       * @pre `imax` and `high` are finite and @f$0 \le high \le imax@f$.
       */
      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::circleState(RealT imax, RealT high)
      {
        const RealT rhs = (imax - high) * (imax + high);
        if (!std::isfinite(rhs) || rhs < ZERO<RealT>)
        {
          return std::numeric_limits<RealT>::quiet_NaN();
        }
        if (rhs == ZERO<RealT>)
        {
          return ZERO<RealT>;
        }

        const RealT ratio = INITIALIZATION_TOLERANCE / rhs;
        return std::sqrt(rhs)
               * std::sqrt(TWO<RealT>
                           / (std::hypot(ratio, TWO<RealT>) + ratio));
      }

      /**
       * @brief Compute off-axis capacity from a continuation state
       *
       * This expression must match the `ilcap` calculation during residual
       * evaluation so initialization lands on the implemented model.
       *
       * @param[in] ilmax Current-circle continuation state on the component base.
       * @return Available off-axis current on the component base.
       */
      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::capacity(RealT ilmax)
      {
        const RealT ilnorm = std::sqrt(ilmax * ilmax + INITIALIZATION_TOLERANCE);
        return (ilmax / ilnorm) * ilmax;
      }

      /**
       * @brief Bisect an initial-limit interval to machine rounding
       *
       * The upper endpoint is returned because initialization needs the first
       * not-below point; solveInitialLimit() separately validates that point as
       * finite and feasible. A finite floating-point interval contains finitely
       * many representable values, so the loop terminates when no interior
       * midpoint remains.
       *
       * @tparam FuncT Monotone predicate type.
       * @param[in] a Lower endpoint, where `below` is true.
       * @param[in] b Upper endpoint, where `below` is false.
       * @param[in] below Predicate returning true only when finite reconstructed
       *                  capacity lies below the requirement.
       * @return The first representable upper-side endpoint.
       * @pre `a` and `b` are finite, `a < b`, and `below` is monotone.
       */
      template <typename scalar_type, typename index_type>
      template <typename FuncT>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::bisect(RealT a, RealT b, FuncT below)
      {
        while (true)
        {
          const RealT mid = std::midpoint(a, b);
          if (mid <= a || b <= mid)
          {
            break;
          }

          if (below(mid))
          {
            a = mid;
          }
          else
          {
            b = mid;
          }
        }

        return b;
      }

      /**
       * @brief Solve the smallest initial total-current limit
       *
       * @param[in] lower Lower bound for the component-base total-current limit.
       * @param[in] high Priority-axis current command on the component base.
       * @param[in] low Required off-axis capacity on the component base.
       * @return The mutually consistent total limit, continuation state, and
       *         off-axis capacity, or `std::nullopt` when no finite feasible
       *         limit is found.
       * @pre The arguments are finite and nonnegative, with `lower >= high` and
       *      `lower >= low`.
       * @warning This function contains conditional branching and may be used
       *          during initialization, but not during residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      auto Reecb<scalar_type, index_type>::solveInitialLimit(
          RealT lower, RealT high, RealT low) -> std::optional<InitialCurrentLimit>
      {
        const RealT ilmax = circleState(lower, high);
        const RealT ilcap = capacity(ilmax);
        if (!std::isfinite(ilcap))
        {
          return std::nullopt;
        }
        if (!(ilcap < low))
        {
          return InitialCurrentLimit{lower, ilmax, ilcap};
        }

        const auto below = [high, low](RealT limit)
        {
          const RealT state = circleState(limit, high);
          const RealT cap   = capacity(state);
          return std::isfinite(cap) && cap < low;
        };

        // Invert ilcap = ilmax^2 / sqrt(ilmax^2 + tolerance), then recover
        // imax from imax^2 = high^2 + ilmax * sqrt(ilmax^2 + tolerance).
        const RealT delta = std::sqrt(INITIALIZATION_TOLERANCE);
        const RealT ilreq = std::sqrt(low)
                            * std::sqrt(HALF<RealT>
                                        * (low + std::hypot(low, TWO<RealT> * delta)));
        const RealT seed = std::hypot(
            high, std::sqrt(ilreq) * std::sqrt(std::hypot(ilreq, delta)));

        const RealT maximum = std::numeric_limits<RealT>::max();
        RealT       a       = lower;
        RealT       b       = seed;
        if (!std::isfinite(b))
        {
          b = maximum;
        }
        else
        {
          b = std::max(b, std::nextafter(a, maximum));
        }
        if (!(a < b))
        {
          return std::nullopt;
        }

        while (below(b))
        {
          a = b;
          if (b >= maximum)
          {
            return std::nullopt;
          }
          if (b > maximum / TWO<RealT>)
          {
            b = maximum;
          }
          else
          {
            b *= TWO<RealT>;
          }
        }

        const RealT result      = bisect(a, b, below);
        const RealT final_state = circleState(result, high);
        const RealT final_cap   = capacity(final_state);
        if (!std::isfinite(result) || !std::isfinite(final_state)
            || !std::isfinite(final_cap) || final_cap < low)
        {
          return std::nullopt;
        }
        return InitialCurrentLimit{result, final_state, final_cap};
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
       * @brief Static method to log time constant warnings
       *
       * @note Used in combination with static std:once_flag and std:call_once,
       *       to reduce the number of times the warning is printed.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "Reecb: any of Trv, Tp, Tiq, or Tpord below "
                       << TIME_CONSTANT_MINIMUM
                       << " s is raised to that floor to keep the controller lags well posed\n";
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
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        va_component_base_ = mva_base_ * static_cast<RealT>(1.0e6);

        if (PfFlag_ && QFlag_)
        {
          Log::warning() << "Reecb: PfFlag and QFlag are both enabled; "
                         << "this is an atypical control configuration\n";
        }

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
      Reecb<scalar_type, index_type>::logOneMinusExp(RealT x)
      {
        static constexpr auto log_two = std::numbers::ln2_v<RealT>;

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      /**
       * @brief Recover the input that produces a requested smooth-clamp output
       *
       * Exact bounds use a finite offset derived from the initialization
       * tolerance; collapsed bounds are reproduced directly.
       *
       * @param[in] output Requested output.
       * @param[in] lower Lower smooth-clamp limit.
       * @param[in] upper Upper smooth-clamp limit.
       * @param[out] input Recovered clamp input on success.
       * @return true when a finite admissible input was recovered.
       * @warning This function contains conditional branching and as such can
       * be used in initialization methods but not in residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      bool Reecb<scalar_type, index_type>::iclamp(RealT output, RealT lower, RealT upper, RealT& input) const
      {
        if (!std::isfinite(output) || !std::isfinite(lower) || !std::isfinite(upper) || lower > upper
            || output < lower - INITIALIZATION_TOLERANCE || output > upper + INITIALIZATION_TOLERANCE)
        {
          return false;
        }

        output = std::clamp(output, lower, upper);
        if (upper == lower)
        {
          input = lower;
          return true;
        }

        const RealT mu     = Math::MU<RealT>;
        const RealT offset = -std::log(std::expm1(mu * HALF<RealT> * INITIALIZATION_TOLERANCE)) / mu;
        if (output == lower)
        {
          input = lower - offset;
          return true;
        }
        if (output == upper)
        {
          input = upper + offset;
          return true;
        }

        const RealT a = mu * (output - lower);
        const RealT b = mu * (upper - output);
        input         = lower + (a + logOneMinusExp(a) - logOneMinusExp(b)) / mu;
        return std::isfinite(input);
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
