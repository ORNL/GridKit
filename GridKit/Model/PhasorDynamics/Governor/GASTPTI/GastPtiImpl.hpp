/**
 * @file GastPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the GASTPTI governor model.
 */

#pragma once

#include <algorithm>
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
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct a GASTPTI governor without parameters
       *
       * The model is sized and every parameter keeps its documented default,
       * which is a valid Normal-mode configuration. No monitor is created, so
       * getMonitor() returns nullptr until the data constructor is used
       * instead.
       */
      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::GastPti()
      {
        size_ = static_cast<IdxT>(GastPtiIdx::MAXIMUM);
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
        size_ = static_cast<IdxT>(GastPtiIdx::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::~GastPti()
      {
      }

      /**
       * @brief Resolve the parameter-derived constants
       *
       * Raises each turbine lag to the well-posedness floor, sizes the
       * component power base from the turbine rating, and encodes the Fixed
       * response mode as the multiplicative mask the residual applies to the
       * three governor states, so the residual keeps a fixed structure for
       * sparse automatic differentiation.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::setDerivedParameters()
      {
        // The lags are raised to the floor in place, so a negative value is
        // rejected here while the value as read is still available. verify()
        // reports the count.
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "GastPti: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(T1_, "T1");
        check_non_negative(T2_, "T2");
        check_non_negative(T3_, "T3");

        if (T1_ < TIME_CONSTANT_MINIMUM || T2_ < TIME_CONSTANT_MINIMUM
            || T3_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "GastPti: T1, T2, and T3 below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the turbine lags well posed\n";
        }

        T1_ = std::max(T1_, TIME_CONSTANT_MINIMUM);
        T2_ = std::max(T2_, TIME_CONSTANT_MINIMUM);
        T3_ = std::max(T3_, TIME_CONSTANT_MINIMUM);

        va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);

        sfix_ = ONE<RealT>;
        if (mode_ == ResponseMode::Fixed)
        {
          sfix_ = ZERO<RealT>;
        }
      }

      /**
       * @brief Convert a system-base power to GASTPTI component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the component base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type GastPti<scalar_type, index_type>::toComponentBase(scalar_type value) const
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
      scalar_type GastPti<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Every parameter is optional and keeps the default documented in the
       * model README when omitted. A non-numeric value, or a response mode
       * that is not an integer, is counted and reported by verify() rather
       * than throwing. Integer JSON values are accepted for real parameters.
       * A Down Only response mode is accepted with a warning and currently
       * simulated as Normal.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
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
        };

        load_real(Params::R, R_, "R");
        load_real(Params::T1, T1_, "T1");
        load_real(Params::T2, T2_, "T2");
        load_real(Params::T3, T3_, "T3");
        load_real(Params::At, At_, "At");
        load_real(Params::Kt, Kt_, "Kt");
        load_real(Params::Vmax, Vmax_, "Vmax");
        load_real(Params::Vmin, Vmin_, "Vmin");
        load_real(Params::Dturb, Dturb_, "Dturb");
        load_real(Params::Trate, Trate_, "Trate");

        if (data.parameters.contains(Params::mode))
        {
          const auto& value = data.parameters.at(Params::mode);
          if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            mode_ = static_cast<ResponseMode>(*index_value);
          }
          else
          {
            Log::error() << "GastPti: parameter 'mode' must be an integer\n";
            ++parameter_error_count_;
          }
        }

        if (mode_ == ResponseMode::DownOnly)
        {
          Log::warning() << "GastPti: mode 2 (Down Only) is not yet supported; "
                            "using mode 0 (Normal)\n";
        }

        setDerivedParameters();
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
       * @brief Bind the monitorable variables to their internal states
       *
       * The mechanical-power output is published on the system base and the
       * remaining outputs on the component base, as documented in the model
       * README.
       */
      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeMonitor()
      {
        using I        = GastPtiIdx;
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::pmech, [this]
                      { return y_.getData()[I::PMECH]; });
        monitor_->set(Variable::xvalve, [this]
                      { return y_.getData()[I::XVALVE]; });
        monitor_->set(Variable::xflow, [this]
                      { return y_.getData()[I::XFLOW]; });
        monitor_->set(Variable::xtemp, [this]
                      { return y_.getData()[I::XTEMP]; });
        monitor_->set(Variable::vload, [this]
                      { return y_.getData()[I::VLOAD]; });
        monitor_->set(Variable::vtemp, [this]
                      { return y_.getData()[I::VTEMP]; });
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return int 0 on success.
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
       * GASTPTI attaches to no bus, so the bus-interface buffer stays empty.
       * Repeated calls reuse the allocated vectors.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::allocate()
      {
        using I = GastPtiIdx;
        using E = GastPtiExt;

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        auto signal_size = E::MAXIMUM;
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<GastPtiInternalVariables::PMECH>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<GastPtiInternalVariables::PMECH>()->set(
              &y[I::PMECH],
              &(this->getVariableIndex(static_cast<IdxT>(I::PMECH))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the GASTPTI configuration
       *
       * Checks parameter-loading errors, static parameter relationships, the
       * mode-dependent valve-limit ordering, and attached external signals.
       * Narrow valve limits that attenuate the smooth anti-windup gate warn
       * without failing verification. Seed feasibility is operating-point
       * dependent and is checked by initialize().
       *
       * @return int Number of configuration errors; zero when valid.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "GastPti: " << message << '\n';
            ret += 1;
          }
        };

        check(R_ > ZERO<RealT>, "R must be positive");
        check(At_ >= ZERO<RealT>, "At must be non-negative");
        check(Kt_ >= ZERO<RealT>, "Kt must be non-negative");

        if (mode_ == ResponseMode::Fixed)
        {
          check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax in Fixed mode");
        }
        else
        {
          check(Vmin_ < Vmax_, "Vmin must be less than Vmax in Normal mode");
        }

        if (mode_ != ResponseMode::Fixed && Vmin_ < Vmax_)
        {
          constexpr RealT minimum_gate{0.99};

          const RealT midpoint = HALF<RealT> * (Vmin_ + Vmax_);
          const RealT maximum_gate =
              Math::indicator(midpoint, ZERO<RealT>, Vmin_, Vmax_);

          if (maximum_gate < minimum_gate)
          {
            Log::warning()
                << "GastPti: narrow valve limits reduce the maximum anti-windup gate to "
                << maximum_gate << "; Normal-mode dynamics may be attenuated\n";
          }
        }

        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(Trate_ > ZERO<RealT>, "Trate must be positive");
        const bool valid_mode = mode_ == ResponseMode::Normal
                                || mode_ == ResponseMode::Fixed
                                || mode_ == ResponseMode::DownOnly;
        check(valid_mode, "mode must be 0 (Normal), 1 (Fixed), or 2 (Down Only)");

        // An attached port must resolve to readable signal storage. The
        // enumerator is a template argument, so each port names itself once.
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

        return ret;
      }

      /**
       * @brief Initialize GASTPTI from the seeded mechanical-power port
       *
       * Reads the assigned system-base `pmech` node and the attached speed
       * input, seats the three turbine states and the LV gate at the seeded
       * fuel flow, solves the load demand that reproduces that flow through
       * the smooth LV gate, and publishes the resolved load reference to an
       * attached `pref` signal. A non-Fixed-mode fuel flow outside the valve
       * limits warns and initializes at the dispatched value. All
       * operating-point checks are completed before model or signal storage
       * is modified.
       *
       * @pre allocate() has completed.
       * @pre The machine model has seeded the assigned `pmech` node.
       *
       * @return int 0 on success; nonzero when the configuration is invalid
       *             or the temperature-gate margin is not positive.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::initialize()
      {
        using I = GastPtiIdx;

        if (verify() > 0)
        {
          Log::error() << "GastPti: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto* y = y_.getData();

        // The assigned pmech node aliases this entry after allocate(). Its
        // system-base seed remains untouched throughout initialization.
        const ScalarT pmech0 = toComponentBase(y[I::PMECH]);

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>();
        }

        const ScalarT xflow0 = pmech0 + Dturb_ * omega0;
        const ScalarT vtemp0 = At_ + Kt_ * (At_ - xflow0);

        if (mode_ != ResponseMode::Fixed
            && (static_cast<RealT>(xflow0) < Vmin_ || static_cast<RealT>(xflow0) > Vmax_))
        {
          Log::warning() << "GastPti: dispatched fuel flow is outside the valve limits "
                            "[Vmin, Vmax]; initializing at the dispatched value\n";
        }

        const RealT margin = static_cast<RealT>(vtemp0 - xflow0);
        if (margin <= ZERO<RealT>)
        {
          Log::error() << "GastPti: initial temperature-gate margin must be positive\n";
          return 1;
        }

        // The LV gate is the smooth minimum, so the load demand must sit
        // where min(VLOAD, VTEMP) reproduces the initial fuel flow exactly.
        const ScalarT vload0 = vtemp0 - Math::iramp(margin);
        const ScalarT pref0  = toSystemBase(vload0 + omega0 / R_);

        y[I::XVALVE] = xflow0;
        y[I::XFLOW]  = xflow0;
        y[I::XTEMP]  = xflow0;
        y[I::VLOAD]  = vload0;
        y[I::VTEMP]  = vtemp0;
        y[I::VLV]    = xflow0;

        pref_set_ = pref0;
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
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::tagDifferentiable()
      {
        using I = GastPtiIdx;

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[I::XVALVE] = true;
        tag_[I::XFLOW]  = true;
        tag_[I::XTEMP]  = true;
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
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Internal residual
       *
       * Evaluates the three turbine states and the four algebraic rows
       * documented in the model README. The body is kept free of branches
       * and loops so that sparse automatic differentiation resolves a fixed
       * structure; the Fixed response mode enters as the multiplicative mask
       * set by setDerivedParameters().
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] wb Bus voltage components; unused, GASTPTI attaches to no bus.
       * @param[in] ws External signal values on system base.
       * @param[out] f Internal residuals.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      GastPti<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        using I = GastPtiIdx;
        using E = GastPtiExt;

        const ScalarT xvalve = y[I::XVALVE];
        const ScalarT xflow  = y[I::XFLOW];
        const ScalarT xtemp  = y[I::XTEMP];
        const ScalarT vload  = y[I::VLOAD];
        const ScalarT vtemp  = y[I::VTEMP];
        const ScalarT vlv    = y[I::VLV];
        const ScalarT pmech  = y[I::PMECH];

        const ScalarT xvalve_dot = yp[I::XVALVE];
        const ScalarT xflow_dot  = yp[I::XFLOW];
        const ScalarT xtemp_dot  = yp[I::XTEMP];

        const ScalarT omega = ws[E::OMEGA];
        const ScalarT pref  = toComponentBase(ws[E::PREF]);

        const ScalarT valve_target = Math::antiwindup(xvalve, vlv - xvalve, Vmin_, Vmax_);

        f[I::XVALVE] = -xvalve_dot + sfix_ * valve_target / T1_;
        f[I::XFLOW]  = -xflow_dot + sfix_ * (-xflow + xvalve) / T2_;
        f[I::XTEMP]  = -xtemp_dot + sfix_ * (-xtemp + xflow) / T3_;
        f[I::VLOAD]  = -omega + R_ * (pref - vload);
        f[I::VTEMP]  = -vtemp + At_ + Kt_ * (At_ - xtemp);
        f[I::VLV]    = -vlv + Math::min(vload, vtemp);
        f[I::PMECH]  = -toComponentBase(pmech) + xflow - Dturb_ * omega;

        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Refreshes the signal interface buffers and evaluates the internal
       * residual. GASTPTI attaches to no bus, so there is no bus interface to
       * refresh. An unattached reference port falls back to the value latched
       * by initialize(); an unattached speed port reads zero deviation.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateResidual()
      {
        using E = GastPtiExt;

        ws_[E::OMEGA] = ZERO<RealT>;
        ws_[E::PREF]  = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<GastPtiExternalVariables::OMEGA>())
        {
          ws_[E::OMEGA] = signals_.template readExternalVariable<GastPtiExternalVariables::OMEGA>();
          ws_indices_[E::OMEGA] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<GastPtiExternalVariables::PREF>())
        {
          ws_[E::PREF] = signals_.template readExternalVariable<GastPtiExternalVariables::PREF>();
          ws_indices_[E::PREF] =
              signals_.template readExternalVariableIndex<GastPtiExternalVariables::PREF>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();
        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
