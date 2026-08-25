/**
 * @file GastPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the GASTPTI governor model.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <mutex>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /// Logger used for GASTPTI diagnostics.
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct a GASTPTI governor without parameters
       *
       * The model is sized and every parameter keeps its documented default.
       * No monitor or mechanical-power output assignment is created, so
       * getMonitor() returns nullptr and verify() rejects the model until the
       * containing system assigns `pmech`.
       */
      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::GastPti()
      {
        size_ = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

      /**
       * @brief Construct a GASTPTI governor from model data
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::GastPti(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(GastPtiInternalVariables::MAXIMUM);
      }

      /**
       * @brief Destroy the governor and its optional variable monitor.
       */
      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::~GastPti()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the mechanical-power output
       *
       * Sizes the state, residual, and signal-interface buffers, seeds the
       * identity index maps, and points the assigned `pmech` node at the
       * internal state it publishes. That node aliases GASTPTI storage from
       * here on, which is how initialize() reads the seed the machine wrote.
       * Repeated allocation by this model reuses its existing link.
       * GASTPTI attaches to no bus, so the bus-interface buffer stays empty.
       *
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::allocate()
      {
        const auto PMECH = static_cast<size_t>(GastPtiInternalVariables::PMECH);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        const auto signal_size = static_cast<size_t>(GastPtiExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (signals_.template isAssigned<GastPtiInternalVariables::PMECH>())
        {
          signals_.template getSignalNode<GastPtiInternalVariables::PMECH>()->set(
              &y[PMECH],
              &(this->getVariableIndex(static_cast<IdxT>(PMECH))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the GASTPTI configuration
       *
       * Checks parameter-loading and time-floor errors, finiteness and static
       * parameter relationships, system/component bases, the required
       * mechanical-power output assignment, attached external signals, and
       * distinct indexed ports.
       *
       * @return int Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::verify() const
      {
        const auto PMECH = static_cast<size_t>(GastPtiInternalVariables::PMECH);

        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "GastPti: " << message << '\n';
            ret += 1;
          }
        };

        check(std::isfinite(R_) && R_ > ZERO<RealT>, "R must be finite and positive");
        check(std::isfinite(At_) && At_ >= ZERO<RealT>,
              "At must be finite and non-negative");
        check(std::isfinite(Kt_) && Kt_ >= ZERO<RealT>,
              "Kt must be finite and non-negative");

        const bool finite_limits = std::isfinite(Vmin_) && std::isfinite(Vmax_);
        check(finite_limits, "Vmin and Vmax must be finite");
        if (finite_limits)
        {
          check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        }

        check(std::isfinite(Dturb_) && Dturb_ >= ZERO<RealT>,
              "Dturb must be finite and non-negative");
        const RealT component_power_base = trate_provided_ ? va_component_base_ : va_system_base_;
        const bool  valid_component_base = std::isfinite(component_power_base)
                                          && component_power_base > ZERO<RealT>;
        const bool valid_system_base =
            std::isfinite(va_system_base_) && va_system_base_ > ZERO<RealT>;
        if (trate_provided_)
        {
          check(std::isfinite(Trate_) && Trate_ > ZERO<RealT>
                    && valid_component_base,
                "Trate must define a finite positive component power base");
        }
        check(valid_system_base, "system power base must be finite and positive");

        if (valid_component_base && valid_system_base)
        {
          const RealT system_to_component = va_system_base_ / component_power_base;
          const RealT component_to_system = component_power_base / va_system_base_;
          check(std::isfinite(system_to_component)
                    && system_to_component > ZERO<RealT>
                    && std::isfinite(component_to_system)
                    && component_to_system > ZERO<RealT>,
                "system/component power-base conversion ratios must be finite and positive");
        }

        check(signals_.template isAssigned<GastPtiInternalVariables::PMECH>(),
              "pmech output must be assigned");

        // An attached port must resolve to writable signal storage.
        auto check_attached_signal =
            [&]<GastPtiExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "GastPti: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<GastPtiExternalVariables::OMEGA>("speed");
        check_attached_signal.template operator()<GastPtiExternalVariables::PREF>("pref");

        const bool omega_linked =
            signals_.template isAttached<GastPtiExternalVariables::OMEGA>()
            && signals_.template isLinked<GastPtiExternalVariables::OMEGA>();
        const bool pref_linked =
            signals_.template isAttached<GastPtiExternalVariables::PREF>()
            && signals_.template isLinked<GastPtiExternalVariables::PREF>();

        if (variable_indices_.size() == static_cast<size_t>(size_))
        {
          const IdxT pmech_index = variable_indices_[PMECH];

          if (omega_linked)
          {
            check(signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>()
                      != pmech_index,
                  "speed and pmech ports must use distinct signals");
          }
          if (pref_linked)
          {
            check(signals_.template readExternalVariableIndex<GastPtiExternalVariables::PREF>()
                      != pmech_index,
                  "pref and pmech ports must use distinct signals");
          }
          if (omega_linked && pref_linked)
          {
            const IdxT omega_index =
                signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>();
            const IdxT pref_index =
                signals_.template readExternalVariableIndex<GastPtiExternalVariables::PREF>();
            if (omega_index != INVALID_INDEX<IdxT>
                || pref_index != INVALID_INDEX<IdxT>)
            {
              check(omega_index != pref_index,
                    "speed and pref ports must use distinct signals");
            }
          }
        }

        return ret;
      }

      /**
       * @brief Initialize GASTPTI from the seeded mechanical-power port
       *
       * Reads the required system-base `pmech` seed and optional speed input,
       * resolves the response limits and valve mask, and seats every
       * algebraic row. For an active response interval, it inverts the smooth
       * low-value selector so its output reproduces the initialized flow; this
       * requires a positive temperature margin. A collapsed response interval
       * holds the valve at its initial position, so its algebraic selector is
       * seated without that active-flow constraint.
       * Every seed, candidate, response bound, and mask is checked before
       * state, response limits, latches, derivatives, or attached signals are
       * modified.
       *
       * @pre allocate() has completed.
       * @pre The machine model has seeded the assigned `pmech` node.
       *
       * @return int 0 on success; nonzero when allocation, configuration,
       *             seed, candidate, response-limit, or active temperature-
       *             gate checks fail.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::initialize()
      {
        const auto XVALVE = static_cast<size_t>(GastPtiInternalVariables::XVALVE);
        const auto XFLOW  = static_cast<size_t>(GastPtiInternalVariables::XFLOW);
        const auto XTEMP  = static_cast<size_t>(GastPtiInternalVariables::XTEMP);
        const auto VLOAD  = static_cast<size_t>(GastPtiInternalVariables::VLOAD);
        const auto VTEMP  = static_cast<size_t>(GastPtiInternalVariables::VTEMP);
        const auto VLV    = static_cast<size_t>(GastPtiInternalVariables::VLV);
        const auto PMECH  = static_cast<size_t>(GastPtiInternalVariables::PMECH);

        if (!allocated_)
        {
          Log::error() << "GastPti: allocate must complete before initialize\n";
          return 1;
        }

        if (verify() > 0)
        {
          Log::error() << "GastPti: cannot initialize with invalid configuration\n";
          return 1;
        }

        if (!trate_provided_)
        {
          va_component_base_ = va_system_base_;
        }

        auto*         y             = y_.getData();
        const ScalarT pmech_system0 = y[PMECH];
        if (!std::isfinite(static_cast<RealT>(pmech_system0)))
        {
          Log::error() << "GastPti: initial pmech seed must be finite\n";
          return 1;
        }

        RealT omega0 = ZERO<RealT>;
        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          omega0 = static_cast<RealT>(
              signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>());
        }
        if (!std::isfinite(omega0))
        {
          Log::error() << "GastPti: initial speed input must be finite\n";
          return 1;
        }

        const ScalarT pmech0       = toComponentBase(pmech_system0);
        const RealT   pmech0_value = static_cast<RealT>(pmech0);
        const RealT   xflow0       = pmech0_value + Dturb_ * omega0;
        const RealT   vtemp0       = At_ + Kt_ * (At_ - xflow0);
        if (!std::isfinite(pmech0_value) || !std::isfinite(xflow0)
            || !std::isfinite(vtemp0))
        {
          Log::error() << "GastPti: initial power, fuel-flow, and temperature candidates must be finite\n";
          return 1;
        }

        const RealT Vmin_response = std::min(Vmin_, xflow0);
        const RealT Vmax_response = std::max(Vmax_, xflow0);

        if (!std::isfinite(Vmin_response) || !std::isfinite(Vmax_response)
            || Vmin_response > Vmax_response)
        {
          Log::error() << "GastPti: initial response limits must be finite and ordered\n";
          return 1;
        }

        if (xflow0 < Vmin_ || xflow0 > Vmax_)
        {
          Log::warning() << "GastPti: initial fuel flow is outside [Vmin, Vmax]; "
                            "response limits are adjusted to include the initialized value\n";
        }

        const bool valve_active = Vmin_response < Vmax_response;
        RealT      vload0       = xflow0;
        RealT      vlv0         = ZERO<RealT>;
        if (valve_active)
        {
          const RealT margin = vtemp0 - xflow0;
          if (!std::isfinite(margin) || margin <= ZERO<RealT>)
          {
            Log::error()
                << "GastPti: active initial temperature-gate margin must be finite and positive\n";
            return 1;
          }

          // Invert the smooth LV gate so its output reproduces the active
          // initial fuel flow. This stable form avoids subtracting two large,
          // nearly equal values when the temperature margin is large.
          vload0 = xflow0 + (margin - iramp(margin));
          vlv0   = xflow0;
        }
        else
        {
          // A fixed valve does not require the selector output to equal xflow.
          // Seat the selector and reference equations at the operating point.
          vlv0 = static_cast<RealT>(Math::min(vload0, vtemp0));
        }

        const RealT pref_component0 = vload0 + omega0 / R_;
        const RealT pref0           = toSystemBase(pref_component0);
        if (!std::isfinite(vload0) || !std::isfinite(vlv0)
            || !std::isfinite(pref_component0) || !std::isfinite(pref0))
        {
          Log::error()
              << "GastPti: initial selector, load-demand, and reference candidates must be finite\n";
          return 1;
        }

        Vmin_response_ = Vmin_response;
        Vmax_response_ = Vmax_response;
        s_valve_       = ZERO<RealT>;
        if (valve_active)
        {
          s_valve_ = ONE<RealT>;
        }

        y[XVALVE] = static_cast<ScalarT>(xflow0);
        y[XFLOW]  = static_cast<ScalarT>(xflow0);
        y[XTEMP]  = static_cast<ScalarT>(xflow0);
        y[VLOAD]  = static_cast<ScalarT>(vload0);
        y[VTEMP]  = static_cast<ScalarT>(vtemp0);
        y[VLV]    = static_cast<ScalarT>(vlv0);

        pref_set_ = static_cast<ScalarT>(pref0);
        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<GastPtiExternalVariables::PREF>(pref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The fuel-valve, fuel-flow, and exhaust-temperature feedback states
       * carry derivatives; every other internal variable is algebraic.
       *
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::tagDifferentiable()
      {
        const auto XVALVE = static_cast<size_t>(GastPtiInternalVariables::XVALVE);
        const auto XFLOW  = static_cast<size_t>(GastPtiInternalVariables::XFLOW);
        const auto XTEMP  = static_cast<size_t>(GastPtiInternalVariables::XTEMP);

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[XVALVE] = true;
        tag_[XFLOW]  = true;
        tag_[XTEMP]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * All GASTPTI variables are per-unit quantities of the same order, so
       * their absolute and relative tolerances use the same value.
       *
       * @param[in] rel_tol Solver relative tolerance.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Evaluate the seven GASTPTI-owned residual rows
       *
       * Refreshes the signal interface buffers and evaluates the internal
       * residual. GASTPTI attaches to no bus, so there is no bus interface to
       * refresh. An unattached reference port falls back to the value latched
       * by initialize(); an unattached speed port reads zero deviation.
       *
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateResidual()
      {
        const auto OMEGA = static_cast<size_t>(GastPtiExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(GastPtiExternalVariables::PREF);

        ws_[OMEGA] = ZERO<RealT>;
        ws_[PREF]  = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<GastPtiExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::PREF>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the GASTPTI signal interface.
       *
       * @return Interface used to assign the mechanical-power output and attach
       *         optional speed and reference inputs.
       */
      template <typename scalar_type, typename index_type>
      auto GastPti<scalar_type, index_type>::getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              GastPtiInternalVariables,
                              GastPtiExternalVariables>&
      {
        return signals_;
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* GastPti<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Internal residual
       *
       * Evaluates the three turbine states and the four algebraic rows
       * documented in the model README. The body is kept free of branches
       * and loops so that sparse automatic differentiation resolves a fixed
       * structure; effective valve limits and the valve mask are resolved
       * before residual evaluation.
       *
       * @param[in] y Internal variables in `GastPtiInternalVariables` order;
       *              each variable uses the base documented by its enum.
       * @param[in] yp Internal derivatives in the same enum order and bases.
       * @param[in] wb Bus voltage components; unused, GASTPTI attaches to no bus.
       * @param[in] ws External signals in `GastPtiExternalVariables` order:
       *               per-unit speed deviation followed by system-base
       *               active-power reference.
       * @param[out] f Model-owned residuals in `GastPtiInternalVariables` order.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline int
      GastPti<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto XVALVE = static_cast<size_t>(GastPtiInternalVariables::XVALVE);
        const auto XFLOW  = static_cast<size_t>(GastPtiInternalVariables::XFLOW);
        const auto XTEMP  = static_cast<size_t>(GastPtiInternalVariables::XTEMP);
        const auto VLOAD  = static_cast<size_t>(GastPtiInternalVariables::VLOAD);
        const auto VTEMP  = static_cast<size_t>(GastPtiInternalVariables::VTEMP);
        const auto VLV    = static_cast<size_t>(GastPtiInternalVariables::VLV);
        const auto PMECH  = static_cast<size_t>(GastPtiInternalVariables::PMECH);

        const auto OMEGA = static_cast<size_t>(GastPtiExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(GastPtiExternalVariables::PREF);

        static_cast<void>(wb);

        const ScalarT xvalve = y[XVALVE];
        const ScalarT xflow  = y[XFLOW];
        const ScalarT xtemp  = y[XTEMP];
        const ScalarT vload  = y[VLOAD];
        const ScalarT vtemp  = y[VTEMP];
        const ScalarT vlv    = y[VLV];
        const ScalarT pmech  = y[PMECH];

        const ScalarT xvalve_dot = yp[XVALVE];
        const ScalarT xflow_dot  = yp[XFLOW];
        const ScalarT xtemp_dot  = yp[XTEMP];

        const ScalarT omega = ws[OMEGA];
        const ScalarT pref  = toComponentBase(ws[PREF]);

        const ScalarT valve_target =
            Math::antiwindup(xvalve, vlv - xvalve, Vmin_response_, Vmax_response_);

        f[XVALVE] = -xvalve_dot + s_valve_ * valve_target / T1_;
        f[XFLOW]  = -xflow_dot + (-xflow + xvalve) / T2_;
        f[XTEMP]  = -xtemp_dot + (-xtemp + xflow) / T3_;
        f[VLOAD]  = -omega + R_ * (pref - vload);
        f[VTEMP]  = -vtemp + At_ + Kt_ * (At_ - xtemp);
        f[VLV]    = -vlv + Math::min(vload, vtemp);
        f[PMECH]  = -toComponentBase(pmech) + xflow - Dturb_ * omega;

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Load one optional real-valued parameter
       *
       * Real and integer serialized values are accepted. Any other stored
       * type records a loading error while preserving the existing default.
       *
       * @param[in] data Model parameter data.
       * @param[in] parameter Parameter key to load.
       * @param[in,out] target Stored parameter value.
       * @param[in] name Serialized parameter name for diagnostics.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::loadRealParameter(
          const ModelDataT& data,
          GastPtiParameters parameter,
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
          Log::error() << "GastPti: parameter '" << name << "' must be numeric\n";
          ++parameter_error_count_;
        }
      }

      /**
       * @brief Validate and floor one turbine time constant
       *
       * Invalid or nonfinite values record a parameter error and are replaced
       * by the floor so later calculations remain well posed. A valid value
       * below the floor is raised and reported through the return value.
       *
       * @param[in,out] value Time constant to validate and floor.
       * @param[in] name Parameter name for diagnostics.
       * @return true when a valid value was raised to the floor.
       */
      template <typename scalar_type, typename index_type>
      bool GastPti<scalar_type, index_type>::floorTimeConstant(
          RealT& value, const char* name)
      {
        if (!std::isfinite(value) || value < ZERO<RealT>)
        {
          Log::error() << "GastPti: " << name
                       << " must be finite and non-negative\n";
          ++parameter_error_count_;
          value = TIME_CONSTANT_MINIMUM;
          return false;
        }

        const bool raised = value < TIME_CONSTANT_MINIMUM;
        value             = std::max(value, TIME_CONSTANT_MINIMUM);
        return raised;
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Every parameter is optional and keeps the default documented in the
       * model README when omitted. A non-numeric value is counted and reported
       * by verify() rather than throwing. Integer JSON values are accepted for
       * real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;
        trate_provided_        = data.parameters.contains(Params::Trate);

        loadRealParameter(data, Params::R, R_, "R");
        loadRealParameter(data, Params::T1, T1_, "T1");
        loadRealParameter(data, Params::T2, T2_, "T2");
        loadRealParameter(data, Params::T3, T3_, "T3");
        loadRealParameter(data, Params::At, At_, "At");
        loadRealParameter(data, Params::Kt, Kt_, "Kt");
        loadRealParameter(data, Params::Vmax, Vmax_, "Vmax");
        loadRealParameter(data, Params::Vmin, Vmin_, "Vmin");
        loadRealParameter(data, Params::Dturb, Dturb_, "Dturb");
        loadRealParameter(data, Params::Trate, Trate_, "Trate");

        setDerivedParameters();
      }

      /**
       * @brief Bind the monitorable variables to their internal states
       *
       * The mechanical-power output is published on the system base and the
       * remaining outputs on the component base, as documented in the model
       * README.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::pmech, [this]
                      { return y_.getData()[static_cast<size_t>(GastPtiInternalVariables::PMECH)]; });
        monitor_->set(Variable::xvalve, [this]
                      { return y_.getData()[static_cast<size_t>(GastPtiInternalVariables::XVALVE)]; });
        monitor_->set(Variable::xflow, [this]
                      { return y_.getData()[static_cast<size_t>(GastPtiInternalVariables::XFLOW)]; });
        monitor_->set(Variable::xtemp, [this]
                      { return y_.getData()[static_cast<size_t>(GastPtiInternalVariables::XTEMP)]; });
        monitor_->set(Variable::vload, [this]
                      { return y_.getData()[static_cast<size_t>(GastPtiInternalVariables::VLOAD)]; });
        monitor_->set(Variable::vtemp, [this]
                      { return y_.getData()[static_cast<size_t>(GastPtiInternalVariables::VTEMP)]; });
      }

      /**
       * @brief Static method to log time constant warnings
       *
       * @note Used in combination with static std:once_flag and std:call_once,
       *       to reduce the number of times the warning is printed. 
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "GastPti: T1, T2, and T3 below "
                       << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor to keep the turbine lags well posed\n";
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Validates and raises each turbine lag in place so every explicit
       * differential row retains a nonzero denominator, sizes the component
       * power base from a supplied turbine rating, and resets the parameter-only
       * response defaults. An omitted rating is resolved from the current system
       * base during initialize(). That method transactionally finalizes the
       * operating-point-dependent response bounds and valve mask.
       * Recording invalid lag inputs before the in-place floor preserves the
       * loading error for verify().
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::setDerivedParameters()
      {
        bool floor_warning = false;

        floor_warning |= floorTimeConstant(T1_, "T1");
        floor_warning |= floorTimeConstant(T2_, "T2");
        floor_warning |= floorTimeConstant(T3_, "T3");

        if (floor_warning)
        {
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        va_component_base_ = ZERO<RealT>;
        if (trate_provided_)
        {
          va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);
        }
        Vmin_response_ = Vmin_;
        Vmax_response_ = Vmax_;

        s_valve_ = ONE<RealT>;
      }

      /**
       * @brief Invert the smooth ramp on its positive range
       *
       * This initialization-only helper is not called by the differentiated
       * residual path.
       *
       * @param[in] value Positive smooth-ramp output.
       * @return Input whose smooth-ramp output equals `value`.
       *
       * @pre `value` is finite and strictly positive.
       */
      template <typename scalar_type, typename index_type>
      auto GastPti<scalar_type, index_type>::iramp(RealT value) -> RealT
      {
        assert(std::isfinite(value) && value > ZERO<RealT>);

        const RealT mu = Math::MU<RealT>;
        return value + std::log(-std::expm1(-mu * value)) / mu;
      }

      /**
       * @brief Convert a system-base power to GASTPTI component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the component base.
       */
      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] inline scalar_type
      GastPti<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * (va_system_base_ / va_component_base_);
      }

      /**
       * @brief Convert a component-base power to the system base
       *
       * @param[in] value Quantity on the component base.
       * @return The same quantity on the system base.
       */
      template <typename scalar_type, typename index_type>
      auto GastPti<scalar_type, index_type>::toSystemBase(RealT value) const -> RealT
      {
        return value * (va_component_base_ / va_system_base_);
      }

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
