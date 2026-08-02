/**
 * @file GastPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the GASTPTI governor model.
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
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
       * @class GastPti
       * @brief Gas turbine governor with speed droop, three turbine lags, and
       *        exhaust-temperature limiting.
       *
       * @tparam scalar_type Plain real or differentiable scalar type.
       * @tparam index_type Integer index type.
       */

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
        size_ = static_cast<IdxT>(index(GastPtiInternalVariables::MAXIMUM));
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
        size_ = static_cast<IdxT>(index(GastPtiInternalVariables::MAXIMUM));
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
        using I = GastPtiInternalVariables;
        using E = GastPtiExternalVariables;

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        auto signal_size = index(E::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<GastPtiInternalVariables::PMECH>())
        {
          auto* pmech = signals_.template getSignalNode<GastPtiInternalVariables::PMECH>();
          auto* y     = y_.getData();
          pmech->set(
              &y[index(I::PMECH)],
              &(this->getVariableIndex(static_cast<IdxT>(index(I::PMECH)))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the GASTPTI configuration
       *
       * Checks parameter-loading and time-floor errors, finiteness and static
       * parameter relationships, the response-mode encoding, system/component
       * bases, the required mechanical-power output assignment, attached
       * external signals, and distinct indexed ports. Operating-point
       * feasibility and response limits that depend on the initial fuel flow
       * are checked by initialize().
       *
       * @return int Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        checkConfiguration(std::isfinite(R_) && R_ > ZERO<RealT>,
                           "R must be finite and positive",
                           ret);
        checkConfiguration(std::isfinite(At_) && At_ >= ZERO<RealT>,
                           "At must be finite and non-negative",
                           ret);
        checkConfiguration(std::isfinite(Kt_) && Kt_ >= ZERO<RealT>,
                           "Kt must be finite and non-negative",
                           ret);

        const bool finite_limits = std::isfinite(Vmin_) && std::isfinite(Vmax_);
        checkConfiguration(finite_limits, "Vmin and Vmax must be finite", ret);
        if (finite_limits)
        {
          checkConfiguration(Vmin_ <= Vmax_,
                             "Vmin must be less than or equal to Vmax",
                             ret);
        }

        checkConfiguration(std::isfinite(Dturb_) && Dturb_ >= ZERO<RealT>,
                           "Dturb must be finite and non-negative",
                           ret);
        const bool valid_component_base = std::isfinite(Trate_)
                                          && Trate_ > ZERO<RealT>
                                          && std::isfinite(va_component_base_)
                                          && va_component_base_ > ZERO<RealT>;
        const bool valid_system_base =
            std::isfinite(va_system_base_) && va_system_base_ > ZERO<RealT>;
        checkConfiguration(valid_component_base,
                           "Trate must define a finite positive component power base",
                           ret);
        checkConfiguration(valid_system_base,
                           "system power base must be finite and positive",
                           ret);

        if (valid_component_base && valid_system_base)
        {
          const RealT system_to_component = va_system_base_ / va_component_base_;
          const RealT component_to_system = va_component_base_ / va_system_base_;
          checkConfiguration(
              std::isfinite(system_to_component)
                  && system_to_component > ZERO<RealT>
                  && std::isfinite(component_to_system)
                  && component_to_system > ZERO<RealT>,
              "system/component power-base conversion ratios must be finite and positive",
              ret);
        }

        const bool valid_mode = mode_ == ResponseMode::Normal
                                || mode_ == ResponseMode::DownOnly
                                || mode_ == ResponseMode::Fixed;
        checkConfiguration(valid_mode,
                           "mode must be 0 (Normal), 1 (Down Only), or 2 (Fixed)",
                           ret);

        checkConfiguration(
            signals_.template isAssigned<GastPtiInternalVariables::PMECH>(),
            "pmech output must be assigned",
            ret);

        const bool omega_attached =
            signals_.template isAttached<GastPtiExternalVariables::OMEGA>();
        const bool pref_attached =
            signals_.template isAttached<GastPtiExternalVariables::PREF>();
        const bool omega_linked =
            omega_attached
            && signals_.template isLinked<GastPtiExternalVariables::OMEGA>();
        const bool pref_linked =
            pref_attached
            && signals_.template isLinked<GastPtiExternalVariables::PREF>();

        if (omega_attached && !omega_linked)
        {
          Log::error() << "GastPti: speed signal attached with no linked source\n";
          ret += 1;
        }
        if (pref_attached && !pref_linked)
        {
          Log::error() << "GastPti: pref signal attached with no linked source\n";
          ret += 1;
        }

        if (variable_indices_.size() == static_cast<size_t>(size_))
        {
          const IdxT pmech_index =
              variable_indices_[index(GastPtiInternalVariables::PMECH)];

          if (omega_linked)
          {
            checkConfiguration(
                signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>()
                    != pmech_index,
                "speed and pmech ports must use distinct signals",
                ret);
          }
          if (pref_linked)
          {
            checkConfiguration(
                signals_.template readExternalVariableIndex<GastPtiExternalVariables::PREF>()
                    != pmech_index,
                "pref and pmech ports must use distinct signals",
                ret);
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
              checkConfiguration(omega_index != pref_index,
                                 "speed and pref ports must use distinct signals",
                                 ret);
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
       * requires a positive temperature margin. Fixed mode and any collapsed
       * response interval hold the valve at its initial position, so their
       * algebraic selector is seated without that active-flow constraint.
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
        using I = GastPtiInternalVariables;

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

        auto*       y             = y_.getData();
        const RealT pmech_system0 = static_cast<RealT>(y[index(I::PMECH)]);
        if (!std::isfinite(pmech_system0))
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

        const RealT pmech0 = pmech_system0 * (va_system_base_ / va_component_base_);
        const RealT xflow0 = pmech0 + Dturb_ * omega0;
        const RealT vtemp0 = At_ + Kt_ * (At_ - xflow0);
        if (!std::isfinite(pmech0) || !std::isfinite(xflow0) || !std::isfinite(vtemp0))
        {
          Log::error() << "GastPti: initial power, fuel-flow, and temperature candidates must be finite\n";
          return 1;
        }

        RealT Vmin_response = Vmin_;
        RealT Vmax_response = Vmax_;
        if (mode_ == ResponseMode::Normal)
        {
          Vmin_response = std::min(Vmin_, xflow0);
          Vmax_response = std::max(Vmax_, xflow0);
        }
        else if (mode_ == ResponseMode::DownOnly)
        {
          Vmin_response = std::min(Vmin_, xflow0);
          Vmax_response = xflow0;
        }
        else
        {
          Vmin_response = xflow0;
          Vmax_response = xflow0;
        }

        if (!std::isfinite(Vmin_response) || !std::isfinite(Vmax_response)
            || Vmin_response > Vmax_response)
        {
          Log::error() << "GastPti: initial response limits must be finite and ordered\n";
          return 1;
        }

        if (mode_ != ResponseMode::Fixed
            && (xflow0 < Vmin_ || xflow0 > Vmax_))
        {
          Log::warning() << "GastPti: initial fuel flow is outside [Vmin, Vmax]; "
                            "response limits are adjusted to include the initialized value\n";
        }

        RealT s_valve = ONE<RealT>;
        if (Vmin_response == Vmax_response)
        {
          s_valve = ZERO<RealT>;
        }

        RealT vload0 = xflow0;
        RealT vlv0   = ZERO<RealT>;
        if (s_valve != ZERO<RealT>)
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
        s_valve_       = s_valve;

        y[index(I::XVALVE)] = static_cast<ScalarT>(xflow0);
        y[index(I::XFLOW)]  = static_cast<ScalarT>(xflow0);
        y[index(I::XTEMP)]  = static_cast<ScalarT>(xflow0);
        y[index(I::VLOAD)]  = static_cast<ScalarT>(vload0);
        y[index(I::VTEMP)]  = static_cast<ScalarT>(vtemp0);
        y[index(I::VLV)]    = static_cast<ScalarT>(vlv0);

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
        using I = GastPtiInternalVariables;

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[index(I::XVALVE)] = true;
        tag_[index(I::XFLOW)]  = true;
        tag_[index(I::XTEMP)]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * All GASTPTI variables are per-unit fuel demands, flows, and powers of
       * the same order, so they share the relative tolerance as their
       * absolute floor.
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
        using E = GastPtiExternalVariables;

        ws_[index(E::OMEGA)] = ZERO<RealT>;
        ws_[index(E::PREF)]  = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          ws_[index(E::OMEGA)] = signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>();
          ws_indices_[index(E::OMEGA)] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          ws_[index(E::PREF)] = signals_.template readExternalVariable<GastPtiExternalVariables::PREF>();
          ws_indices_[index(E::PREF)] =
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
       * structure; response-limit modes enter through valve limits and the
       * valve mask resolved before residual evaluation.
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
        using I = GastPtiInternalVariables;
        using E = GastPtiExternalVariables;

        static_cast<void>(wb);

        const ScalarT xvalve = y[index(I::XVALVE)];
        const ScalarT xflow  = y[index(I::XFLOW)];
        const ScalarT xtemp  = y[index(I::XTEMP)];
        const ScalarT vload  = y[index(I::VLOAD)];
        const ScalarT vtemp  = y[index(I::VTEMP)];
        const ScalarT vlv    = y[index(I::VLV)];
        const ScalarT pmech  = y[index(I::PMECH)];

        const ScalarT xvalve_dot = yp[index(I::XVALVE)];
        const ScalarT xflow_dot  = yp[index(I::XFLOW)];
        const ScalarT xtemp_dot  = yp[index(I::XTEMP)];

        const ScalarT omega = ws[index(E::OMEGA)];
        const ScalarT pref  = toComponentBase(ws[index(E::PREF)]);

        const ScalarT valve_target =
            Math::antiwindup(xvalve, vlv - xvalve, Vmin_response_, Vmax_response_);

        f[index(I::XVALVE)] = -xvalve_dot + s_valve_ * valve_target / T1_;
        f[index(I::XFLOW)]  = -xflow_dot + (-xflow + xvalve) / T2_;
        f[index(I::XTEMP)]  = -xtemp_dot + (-xtemp + xflow) / T3_;
        f[index(I::VLOAD)]  = -omega + R_ * (pref - vload);
        f[index(I::VTEMP)]  = -vtemp + At_ + Kt_ * (At_ - xtemp);
        f[index(I::VLV)]    = -vlv + Math::min(vload, vtemp);
        f[index(I::PMECH)]  = -toComponentBase(pmech) + xflow - Dturb_ * omega;

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Record one failed configuration condition
       *
       * @param[in] condition Required condition.
       * @param[in] message Error message when `condition` is false.
       * @param[in,out] errors Accumulated configuration-error count.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::checkConfiguration(
          bool condition, const char* message, int& errors)
      {
        if (!condition)
        {
          Log::error() << "GastPti: " << message << '\n';
          errors += 1;
        }
      }

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
       * model README when omitted. A non-numeric value, or a response mode
       * that is not an integer, is counted and reported by verify() rather
       * than throwing. Integer JSON values are accepted for real parameters.
       * The serialized response-mode values are decoded explicitly as
       * 0 Normal, 1 Down Only, and 2 Fixed.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

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

        if (data.parameters.contains(Params::mode))
        {
          const auto& value = data.parameters.at(Params::mode);
          if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            switch (*index_value)
            {
            case 0:
              mode_ = ResponseMode::Normal;
              break;
            case 1:
              mode_ = ResponseMode::DownOnly;
              break;
            case 2:
              mode_ = ResponseMode::Fixed;
              break;
            default:
              Log::error() << "GastPti: parameter 'mode' must be 0 (Normal), "
                              "1 (Down Only), or 2 (Fixed)\n";
              ++parameter_error_count_;
              break;
            }
          }
          else
          {
            Log::error() << "GastPti: parameter 'mode' must be an integer\n";
            ++parameter_error_count_;
          }
        }

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
        using I        = GastPtiInternalVariables;
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::pmech, [this]
                      { return y_.getData()[index(I::PMECH)]; });
        monitor_->set(Variable::xvalve, [this]
                      { return y_.getData()[index(I::XVALVE)]; });
        monitor_->set(Variable::xflow, [this]
                      { return y_.getData()[index(I::XFLOW)]; });
        monitor_->set(Variable::xtemp, [this]
                      { return y_.getData()[index(I::XTEMP)]; });
        monitor_->set(Variable::vload, [this]
                      { return y_.getData()[index(I::VLOAD)]; });
        monitor_->set(Variable::vtemp, [this]
                      { return y_.getData()[index(I::VTEMP)]; });
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Validates and raises each turbine lag in place so every explicit
       * differential row retains a nonzero denominator, sizes the component
       * power base from the turbine rating, and resets the parameter-only
       * response defaults. initialize() transactionally finalizes the
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
          Log::warning() << "GastPti: T1, T2, and T3 below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the turbine lags well posed\n";
        }

        va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);
        Vmin_response_     = Vmin_;
        Vmax_response_     = Vmax_;

        s_valve_ = ONE<RealT>;
        if (mode_ == ResponseMode::Fixed)
        {
          s_valve_ = ZERO<RealT>;
        }
      }

      /**
       * @brief Invert the smooth ramp on its positive range
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
