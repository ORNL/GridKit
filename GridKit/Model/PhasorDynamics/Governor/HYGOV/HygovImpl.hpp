/**
 * @file HygovImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the HYGOV governor model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct a HYGOV governor without parameters
       *
       * The model is sized with the documented parameter defaults but without
       * a monitor or assigned mechanical-power output.
       */
      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov()
      {
        size_ = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
      }

      /**
       * @brief Construct a HYGOV governor from model data
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::~Hygov()
      {
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the mechanical-power output
       *
       * Sizes the state, residual, and signal-interface buffers, initializes
       * the identity index maps, and points the assigned `pmech` node at the
       * internal state it publishes. That node aliases HYGOV storage from
       * here on, which is how initialize() reads the machine value.
       * HYGOV attaches to no bus, so the bus-interface buffer stays empty.
       * Repeated calls reuse the allocated vectors.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::allocate()
      {
        const auto PMECH = static_cast<size_t>(HygovInternalVariables::PMECH);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        const auto signal_size = static_cast<size_t>(HygovExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<HygovInternalVariables::PMECH>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<HygovInternalVariables::PMECH>()->set(
              &y[PMECH],
              &(this->getVariableIndex(static_cast<IdxT>(PMECH))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the HYGOV configuration
       *
       * Checks parameter-loading errors, static parameter relationships, the
       * power-base domain, gate-curve shape, gate-limit domain, the assigned
       * mechanical-power output, and attached external signals. Mechanical-
       * power feasibility is operating-point dependent and is checked by
       * initialize().
       *
       * @return int Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Hygov: " << message << '\n';
            ret += 1;
          }
        };

        RealT      component_power_base = va_component_base_;
        const bool component_base_is_omitted =
            !(component_power_base > ZERO<RealT>);
        if (component_base_is_omitted)
        {
          component_power_base = va_system_base_;
        }

        const bool valid_component_base = std::isfinite(component_power_base)
                                          && component_power_base > ZERO<RealT>;
        const bool valid_system_base = std::isfinite(va_system_base_)
                                       && va_system_base_ > ZERO<RealT>;
        check(valid_component_base,
              "component power base must be finite and positive");
        check(valid_system_base,
              "system power base must be finite and positive");
        if (valid_component_base && valid_system_base)
        {
          const RealT system_to_component = va_system_base_ / component_power_base;
          const RealT component_to_system = component_power_base / va_system_base_;
          const bool  valid_base_ratios   = std::isfinite(system_to_component)
                                         && system_to_component > ZERO<RealT>
                                         && std::isfinite(component_to_system)
                                         && component_to_system > ZERO<RealT>;
          check(valid_base_ratios,
                "system/component power-base conversion ratios must be finite and positive");
        }

        check(Rtemp_ < ZERO<RealT> || Rtemp_ > ZERO<RealT>, "Rtemp must be nonzero");
        check(Tn_ >= ZERO<RealT>, "Tn must be non-negative");
        check(Velm_ >= ZERO<RealT>, "Velm must be non-negative");
        check(Gmin_ < Gmax_, "Gmin must be less than Gmax");
        check(At_ > ZERO<RealT>, "At must be positive");
        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(db1_ >= ZERO<RealT>, "db1 must be non-negative");
        check(Hdam_ > ZERO<RealT>, "Hdam must be positive");

        bool curve_shape_is_valid = true;
        for (size_t i = 1; i < Gv_.size(); ++i)
        {
          const bool gate_points_increase  = Gv_[i - 1] < Gv_[i];
          const bool power_points_increase = Pgv_[i - 1] <= Pgv_[i];

          check(gate_points_increase, "Gv points must be strictly increasing");
          check(power_points_increase, "Pgv points must be non-decreasing");

          if (!gate_points_increase || !power_points_increase)
          {
            curve_shape_is_valid = false;
          }
        }
        const bool minimum_gate_is_valid = Gv_[0] <= Gmin_;
        const bool maximum_gate_is_valid = Gmax_ <= Gv_[5];
        check(minimum_gate_is_valid, "Gmin must be at or above the first Gv point");
        check(maximum_gate_is_valid, "Gmax must be at or below the last Gv point");

        const bool can_check_power_range = curve_shape_is_valid
                                           && Gmin_ < Gmax_
                                           && minimum_gate_is_valid
                                           && maximum_gate_is_valid
                                           && At_ > ZERO<RealT>
                                           && Hdam_ > ZERO<RealT>;
        if (can_check_power_range)
        {
          // A rise no wider than the tolerance that pins a seed to a range
          // edge leaves the gate undetermined by the mechanical power.
          const RealT minimum_power      = initialMechanicalPower(Gmin_);
          const RealT maximum_power      = initialMechanicalPower(Gmax_);
          const RealT power_range        = maximum_power - minimum_power;
          const bool  finite_power_range = std::isfinite(minimum_power)
                                          && std::isfinite(maximum_power)
                                          && std::isfinite(power_range);
          check(finite_power_range,
                "mechanical-power range must be finite");
          if (finite_power_range)
          {
            check(power_range > INITIALIZATION_TOLERANCE,
                  "mechanical power must rise across [Gmin, Gmax]");
          }
        }

        check(signals_.template isAssigned<HygovInternalVariables::PMECH>(),
              "pmech output signal must be assigned");

        // An attached port must resolve to readable signal storage. The
        // enumerator is a template argument, so each port names itself once.
        auto check_attached_signal =
            [&]<HygovExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Hygov: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<HygovExternalVariables::OMEGA>("speed");
        check_attached_signal.template operator()<HygovExternalVariables::PREF>("pref");
        check_attached_signal.template operator()<HygovExternalVariables::PAUX>("paux");

        return ret;
      }

      /**
       * @brief Initialize HYGOV from the mechanical-power port
       *
       * Reads the assigned system-base `pmech` node and the attached speed
       * and auxiliary-power inputs, solves the component-base steady state
       * that preserves the given value, and publishes the resolved load reference
       * to an attached `pref` signal.
       *
       * @pre allocate() has completed.
       * @pre The machine model has initialized the assigned `pmech` node.
       *
       * @post On success the state zeros every residual row at machine
       *       rounding; a value clipped to the achievable-power range edge
       *       leaves a mechanical-power residual up to the initialization
       *       tolerance.
       * @post On failure no state or signal storage has changed.
       *
       * @return int 0 on success; nonzero when the configuration or initial
       *             values are invalid, the initial speed deviation is
       *             nonzero, or no gate inside [Gmin, Gmax] reproduces the
       *             given power.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::initialize()
      {
        const auto XN      = static_cast<size_t>(HygovInternalVariables::XN);
        const auto XF      = static_cast<size_t>(HygovInternalVariables::XF);
        const auto C       = static_cast<size_t>(HygovInternalVariables::C);
        const auto G       = static_cast<size_t>(HygovInternalVariables::G);
        const auto Q       = static_cast<size_t>(HygovInternalVariables::Q);
        const auto OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto FC      = static_cast<size_t>(HygovInternalVariables::FC);
        const auto RC      = static_cast<size_t>(HygovInternalVariables::RC);
        const auto PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        const auto H       = static_cast<size_t>(HygovInternalVariables::H);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);

        bool ret = verify() == 0;
        if (!ret)
        {
          Log::error() << "Hygov: cannot initialize with invalid configuration\n";
          return 1;
        }

        ret = va_component_base_ > ZERO<RealT>;
        if (!ret)
        {
          va_component_base_ = va_system_base_;
        }

        auto* y = y_.getData();

        // The assigned pmech node aliases this entry after allocate(). Its
        // system-base value remains untouched throughout initialization.
        const ScalarT pmech0_system = y[PMECH];

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<HygovExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<HygovExternalVariables::OMEGA>();
        }

        ScalarT paux0_system{ZERO<RealT>};
        if (signals_.template isAttached<HygovExternalVariables::PAUX>())
        {
          paux0_system = signals_.template readExternalVariable<HygovExternalVariables::PAUX>();
        }

        auto is_finite = [](ScalarT value)
        {
          return std::isfinite(static_cast<RealT>(value));
        };
        ret = is_finite(pmech0_system)
              && is_finite(omega0)
              && is_finite(paux0_system);
        if (!ret)
        {
          Log::error() << "Hygov: initial pmech, speed, and paux values must be finite\n";
          return 1;
        }

        const ScalarT pmech0 = toComponentBase(pmech0_system);
        const ScalarT paux0  = toComponentBase(paux0_system);
        ret                  = is_finite(pmech0)
              && is_finite(paux0);
        if (!ret)
        {
          Log::error() << "Hygov: initial power-base conversions must be finite\n";
          return 1;
        }

        // Synchronous machines provide an exactly zero speed deviation. A
        // moving machine would need a multi-root gate search, which this
        // model does not support. The speed was verified finite above, so
        // the two-sided sign test is exact.
        const RealT speed0 = static_cast<RealT>(omega0);

        ret = !(speed0 < ZERO<RealT>) && !(speed0 > ZERO<RealT>);
        if (!ret)
        {
          Log::error() << "Hygov: initialization requires zero speed deviation\n";
          return 1;
        }

        const RealT gate0 = solveInitialGate(static_cast<RealT>(pmech0));
        ret               = std::isfinite(gate0);
        if (!ret)
        {
          Log::error()
              << "Hygov: no gate inside [Gmin, Gmax] reproduces the given mechanical power\n";
          return 1;
        }

        const ScalarT h0       = static_cast<ScalarT>(Hdam_);
        const ScalarT pgv0     = gatePower(static_cast<ScalarT>(gate0));
        const ScalarT q0       = std::sqrt(Hdam_) * pgv0;
        const ScalarT omegadb0 = Math::deadband1(omega0, -db1_, db1_);
        const ScalarT xn0      = omegadb0;
        const ScalarT yomega0  = xn0 + leadlag_gain_ * (omegadb0 - xn0);
        const ScalarT pref0    = toSystemBase(yomega0 + Rperm_ * gate0 - paux0);

        ret = is_finite(h0)
              && is_finite(pgv0)
              && is_finite(q0)
              && is_finite(omegadb0)
              && is_finite(xn0)
              && is_finite(yomega0)
              && is_finite(pref0);
        if (!ret)
        {
          Log::error() << "Hygov: initialization produced a nonfinite value\n";
          return 1;
        }

        y[XN]      = xn0;
        y[XF]      = ZERO<RealT>;
        y[C]       = gate0;
        y[G]       = gate0;
        y[Q]       = q0;
        y[OMEGADB] = omegadb0;
        y[EF]      = ZERO<RealT>;
        y[FC]      = ZERO<RealT>;
        y[RC]      = ZERO<RealT>;
        y[PGV]     = pgv0;
        y[H]       = h0;

        pref_set_ = pref0;
        paux_set_ = paux0_system;

        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<HygovExternalVariables::PREF>(pref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Identify the differential variables
       *
       * The speed lead-lag state, the governor error filter, the desired
       * gate, the gate servo, and the turbine flow carry derivatives; every
       * other internal variable is algebraic.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::tagDifferentiable()
      {
        const auto XN = static_cast<size_t>(HygovInternalVariables::XN);
        const auto XF = static_cast<size_t>(HygovInternalVariables::XF);
        const auto C  = static_cast<size_t>(HygovInternalVariables::C);
        const auto G  = static_cast<size_t>(HygovInternalVariables::G);
        const auto Q  = static_cast<size_t>(HygovInternalVariables::Q);

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[XN] = true;
        tag_[XF] = true;
        tag_[C]  = true;
        tag_[G]  = true;
        tag_[Q]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * All HYGOV variables are per-unit speeds, gates, flows, heads, and
       * powers of the same order, so their absolute and relative tolerance
       * have the same value.
       *
       * @param[in] rel_tol Solver relative tolerance.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Refreshes the signal interface buffers and evaluates the internal
       * residual. HYGOV attaches to no bus, so there is no bus interface to
       * refresh. An unattached reference or auxiliary port falls back to the
       * value latched by initialize(); an unattached speed port reads zero
       * deviation.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateResidual()
      {
        const auto OMEGA = static_cast<size_t>(HygovExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(HygovExternalVariables::PREF);
        const auto PAUX  = static_cast<size_t>(HygovExternalVariables::PAUX);

        ws_[OMEGA] = ZERO<RealT>;
        ws_[PREF]  = pref_set_;
        ws_[PAUX]  = paux_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<HygovExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<HygovExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<HygovExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::PREF>();
        }
        if (signals_.template isAttached<HygovExternalVariables::PAUX>())
        {
          ws_[PAUX] = signals_.template readExternalVariable<HygovExternalVariables::PAUX>();
          ws_indices_[PAUX] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::PAUX>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();
        return 0;
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Hygov<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Internal residual
       *
       * Evaluates the five governor states and the seven algebraic rows
       * documented in the model README. The body is kept free of branches
       * and loops so that sparse automatic differentiation resolves a fixed
       * structure; the gate curve enters as a fixed sum of smooth linear
       * segments.
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] wb Bus voltage components; unused, HYGOV attaches to no bus.
       * @param[in] ws External signal values on system base.
       * @param[out] f Internal residuals.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Hygov<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto XN      = static_cast<size_t>(HygovInternalVariables::XN);
        const auto XF      = static_cast<size_t>(HygovInternalVariables::XF);
        const auto C       = static_cast<size_t>(HygovInternalVariables::C);
        const auto G       = static_cast<size_t>(HygovInternalVariables::G);
        const auto Q       = static_cast<size_t>(HygovInternalVariables::Q);
        const auto OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto FC      = static_cast<size_t>(HygovInternalVariables::FC);
        const auto RC      = static_cast<size_t>(HygovInternalVariables::RC);
        const auto PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        const auto H       = static_cast<size_t>(HygovInternalVariables::H);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);

        const auto OMEGA = static_cast<size_t>(HygovExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(HygovExternalVariables::PREF);
        const auto PAUX  = static_cast<size_t>(HygovExternalVariables::PAUX);

        const ScalarT xn      = y[XN];
        const ScalarT xf      = y[XF];
        const ScalarT c       = y[C];
        const ScalarT g       = y[G];
        const ScalarT q       = y[Q];
        const ScalarT omegadb = y[OMEGADB];
        const ScalarT ef      = y[EF];
        const ScalarT fc      = y[FC];
        const ScalarT rc      = y[RC];
        const ScalarT pgv     = y[PGV];
        const ScalarT head    = y[H];
        const ScalarT pmech   = y[PMECH];

        const ScalarT xn_dot = yp[XN];
        const ScalarT xf_dot = yp[XF];
        const ScalarT c_dot  = yp[C];
        const ScalarT g_dot  = yp[G];
        const ScalarT q_dot  = yp[Q];

        const ScalarT omega = ws[OMEGA];
        const ScalarT pref  = ws[PREF];
        const ScalarT paux  = ws[PAUX];

        const ScalarT yomega = xn + leadlag_gain_ * (omegadb - xn);

        f[XN]      = -xn_dot + (omegadb - xn) / Tnp_;
        f[XF]      = -xf_dot + (ef - xf) / Tf_;
        f[C]       = -c_dot + Math::antiwindup(c, rc, Gmin_, Gmax_);
        f[G]       = -g_dot + (c - g) / Tg_;
        f[Q]       = -q_dot + (Hdam_ - head) / Tw_;
        f[OMEGADB] = -omegadb + Math::deadband1(omega, -db1_, db1_);
        f[EF]      = -ef + toComponentBase(pref + paux) - yomega - Rperm_ * c;
        f[FC]      = -Rtemp_ * fc + xf / Tr_ + (ef - xf) / Tf_;
        f[RC]      = -rc + Math::clamp(fc, -Velm_, Velm_);
        f[PGV]     = -pgv + gatePower(g);
        f[H]       = -q * q + head * pgv * pgv;
        f[PMECH]   = -toComponentBase(pmech) + At_ * head * (q - Qnl_) - Dturb_ * omega * g;

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Read the parameters out of the model data
       *
       * Every omitted parameter keeps the default documented in the model
       * README. A non-numeric or nonfinite value is counted and reported by
       * verify() rather than throwing. Integer JSON values are accepted for
       * real parameters. All-zero `Gv` and `Pgv` source points select the
       * identity gate curve.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_real = [&](auto key, RealT& target, const char* name) -> bool
        {
          if (!data.parameters.contains(key))
          {
            return false;
          }

          const auto& value = data.parameters.at(key);
          RealT       parsed_value{};
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            parsed_value = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            parsed_value = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "Hygov: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
            return false;
          }

          const bool ret = std::isfinite(parsed_value);
          if (!ret)
          {
            Log::error() << "Hygov: parameter '" << name << "' must be finite\n";
            ++parameter_error_count_;
            return false;
          }

          target = parsed_value;
          return true;
        };

        bool ret = load_real(Params::Trate, va_component_base_, "Trate");
        if (ret)
        {
          ret = va_component_base_ > ZERO<RealT>;
          if (!ret)
          {
            Log::error() << "Hygov: Trate must be positive when provided\n";
            ++parameter_error_count_;
          }
          va_component_base_ *= static_cast<RealT>(1.0e6);
        }
        load_real(Params::Rperm, Rperm_, "Rperm");
        load_real(Params::Rtemp, Rtemp_, "Rtemp");
        load_real(Params::Tr, Tr_, "Tr");
        load_real(Params::Tf, Tf_, "Tf");
        load_real(Params::Tg, Tg_, "Tg");
        load_real(Params::Velm, Velm_, "Velm");
        load_real(Params::Gmax, Gmax_, "Gmax");
        load_real(Params::Gmin, Gmin_, "Gmin");
        load_real(Params::Tw, Tw_, "Tw");
        load_real(Params::At, At_, "At");
        load_real(Params::Dturb, Dturb_, "Dturb");
        load_real(Params::Qnl, Qnl_, "Qnl");
        load_real(Params::Tn, Tn_, "Tn");
        load_real(Params::Tnp, Tnp_, "Tnp");
        load_real(Params::db1, db1_, "db1");
        load_real(Params::db2, db2_, "db2");
        load_real(Params::Hdam, Hdam_, "Hdam");
        load_real(Params::Gv0, Gv_[0], "Gv0");
        load_real(Params::Gv1, Gv_[1], "Gv1");
        load_real(Params::Gv2, Gv_[2], "Gv2");
        load_real(Params::Gv3, Gv_[3], "Gv3");
        load_real(Params::Gv4, Gv_[4], "Gv4");
        load_real(Params::Gv5, Gv_[5], "Gv5");
        load_real(Params::Pgv0, Pgv_[0], "Pgv0");
        load_real(Params::Pgv1, Pgv_[1], "Pgv1");
        load_real(Params::Pgv2, Pgv_[2], "Pgv2");
        load_real(Params::Pgv3, Pgv_[3], "Pgv3");
        load_real(Params::Pgv4, Pgv_[4], "Pgv4");
        load_real(Params::Pgv5, Pgv_[5], "Pgv5");

        auto deviates_from_zero = [](RealT value)
        { return value < ZERO<RealT> || value > ZERO<RealT>; };

        const bool curve_supplied =
            std::any_of(Gv_.begin(), Gv_.end(), deviates_from_zero)
            || std::any_of(Pgv_.begin(), Pgv_.end(), deviates_from_zero);
        if (!curve_supplied)
        {
          Gv_  = {ZERO<RealT>,
                  static_cast<RealT>(0.2),
                  static_cast<RealT>(0.4),
                  static_cast<RealT>(0.6),
                  static_cast<RealT>(0.8),
                  ONE<RealT>};
          Pgv_ = Gv_;
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
      void Hygov<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::pmech, [this]
                      { return y_.getData()[static_cast<size_t>(HygovInternalVariables::PMECH)]; });
        monitor_->set(Variable::filter, [this]
                      { return y_.getData()[static_cast<size_t>(HygovInternalVariables::XF)]; });
        monitor_->set(Variable::desiredgate, [this]
                      { return y_.getData()[static_cast<size_t>(HygovInternalVariables::C)]; });
        monitor_->set(Variable::gate, [this]
                      { return y_.getData()[static_cast<size_t>(HygovInternalVariables::G)]; });
        monitor_->set(Variable::flow, [this]
                      { return y_.getData()[static_cast<size_t>(HygovInternalVariables::Q)]; });
        monitor_->set(Variable::head, [this]
                      { return y_.getData()[static_cast<size_t>(HygovInternalVariables::H)]; });
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Floors each governor time constant so the residual equations retain
       * Hessenberg form, then derives the speed lead-lag gain.
       */
      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::setDerivedParameters()
      {
        // The lags are raised to the floor in place, so a negative value is
        // rejected here while the value as read is still available. verify()
        // reports the count.
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Hygov: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(Tr_, "Tr");
        check_non_negative(Tf_, "Tf");
        check_non_negative(Tg_, "Tg");
        check_non_negative(Tw_, "Tw");
        check_non_negative(Tnp_, "Tnp");

        if (Tr_ < TIME_CONSTANT_MINIMUM || Tf_ < TIME_CONSTANT_MINIMUM
            || Tg_ < TIME_CONSTANT_MINIMUM || Tw_ < TIME_CONSTANT_MINIMUM
            || Tnp_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Hygov: Tr, Tf, Tg, Tw, and Tnp below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to preserve Hessenberg form\n";
        }

        // HYGOV residuals solve explicitly for the state derivatives to preserve
        // Hessenberg form. A zero time constant would instead require an implicit
        // residual formulation, so enforce a strictly positive lower bound.
        Tr_  = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Tf_  = std::max(Tf_, TIME_CONSTANT_MINIMUM);
        Tg_  = std::max(Tg_, TIME_CONSTANT_MINIMUM);
        Tw_  = std::max(Tw_, TIME_CONSTANT_MINIMUM);
        Tnp_ = std::max(Tnp_, TIME_CONSTANT_MINIMUM);

        leadlag_gain_ = Tn_ / Tnp_;
      }

      /**
       * @brief Evaluate the nonlinear gate-to-power curve
       *
       * Sums the five smooth CommonMath linear segments spanned by the
       * `Gv`/`Pgv` points, so the same fixed expression serves the residual
       * and both scalar instantiations.
       *
       * @param[in] gate Gate position.
       * @return Turbine power at nominal head.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline scalar_type
      Hygov<scalar_type, index_type>::gatePower(scalar_type gate) const
      {
        ScalarT retval = Pgv_[0]
                         + Math::linseg(gate, Gv_[0], Gv_[1], Pgv_[1] - Pgv_[0])
                         + Math::linseg(gate, Gv_[1], Gv_[2], Pgv_[2] - Pgv_[1])
                         + Math::linseg(gate, Gv_[2], Gv_[3], Pgv_[3] - Pgv_[2])
                         + Math::linseg(gate, Gv_[3], Gv_[4], Pgv_[4] - Pgv_[3])
                         + Math::linseg(gate, Gv_[4], Gv_[5], Pgv_[5] - Pgv_[4]);

        return retval;
      }

      /**
       * @brief Steady component-base mechanical power at a gate position
       *
       * At the steady state the head equals the dam head and the flow
       * follows the gate curve, so the PGV, H, and PMECH rows collapse to
       * @f[
       *   P_{\mathrm{m}}(g)
       *     = A_t H_{\mathrm{dam}}
       *       \left(\sqrt{H_{\mathrm{dam}}}\,N_{\mathrm{GV}}(g)
       *             - q_{\mathrm{NL}}\right).
       * @f]
       * The expression is composed exactly as those rows compose it, so a
       * gate solved against it zeros the implemented residual at machine
       * rounding.
       *
       * @param[in] gate Gate position.
       * @return Steady mechanical power on the component base.
       */
      template <typename scalar_type, typename index_type>
      typename Hygov<scalar_type, index_type>::RealT
      Hygov<scalar_type, index_type>::initialMechanicalPower(RealT gate) const
      {
        const RealT pgv = static_cast<RealT>(gatePower(static_cast<ScalarT>(gate)));
        const RealT q   = std::sqrt(Hdam_) * pgv;
        return At_ * Hdam_ * (q - Qnl_);
      }

      /**
       * @brief Solve the steady gate position for a given mechanical power
       *
       * Initialization requires a zero speed deviation and verify() requires
       * the steady power to rise across [Gmin, Gmax], so the endpoint
       * residuals decide feasibility and bisection converges to a root of
       * the nondecreasing steady-power curve.
       *
       * @pre verify() reports no errors.
       *
       * @param[in] pmech Mechanical power on the component base.
       * @return The gate position, or a quiet NaN when no gate inside
       *         [Gmin, Gmax] reproduces the value within the initialization
       *         tolerance.
       *
       * @warning This function contains conditional branching and may be used
       *          during initialization, but not during residual evaluation.
       */
      template <typename scalar_type, typename index_type>
      typename Hygov<scalar_type, index_type>::RealT
      Hygov<scalar_type, index_type>::solveInitialGate(RealT pmech) const
      {
        // A nonfinite seed is unreproducible by any gate and would otherwise
        // slip through the sign tests below.
        const bool ret = std::isfinite(pmech);
        if (!ret)
        {
          return std::numeric_limits<RealT>::quiet_NaN();
        }

        RealT a  = Gmin_;
        RealT b  = Gmax_;
        RealT fa = initialMechanicalPower(a) - pmech;
        RealT fb = initialMechanicalPower(b) - pmech;

        // A value just outside the achievable range pins to the gate limit
        // when it is within the initialization tolerance of the range edge.
        if (fa > ZERO<RealT>)
        {
          if (fa <= INITIALIZATION_TOLERANCE)
          {
            return a;
          }
          return std::numeric_limits<RealT>::quiet_NaN();
        }
        if (fb < ZERO<RealT>)
        {
          if (-fb <= INITIALIZATION_TOLERANCE)
          {
            return b;
          }
          return std::numeric_limits<RealT>::quiet_NaN();
        }

        // Bisect until no representable midpoint remains, then keep the
        // endpoint with the smaller residual. Termination is guaranteed:
        // std::midpoint cannot overflow and lands inside [a, b], so every
        // accepted step strictly shrinks the finite set of representable
        // values between the endpoints.
        while (true)
        {
          const RealT mid = std::midpoint(a, b);
          if (mid <= a || b <= mid)
          {
            break;
          }
          const RealT fmid = initialMechanicalPower(mid) - pmech;
          if (fmid <= ZERO<RealT>)
          {
            a  = mid;
            fa = fmid;
          }
          else
          {
            b  = mid;
            fb = fmid;
          }
        }
        if (std::abs(fa) <= std::abs(fb))
        {
          return a;
        }
        return b;
      }

      /**
       * @brief Convert a system-base power to HYGOV component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the component base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_component_base_;
      }

      /**
       * @brief Convert a component-base power to the system base
       *
       * @param[in] value Quantity on the component base.
       * @return The same quantity on the system base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
