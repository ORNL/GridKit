/**
 * @file ReecbImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REECB electrical-control model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief Construct a REECB controller without parameters
       *
       * The model is sized but left unconfigured. Every parameter keeps its
       * documented default, the required power base is absent, and no monitor
       * is created, so verify() reports configuration errors until the data
       * constructor is used instead.
       *
       * @param[in] bus Terminal bus the controller measures.
       */
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(ReecbIdx::MAXIMUM);
      }

      /**
       * @brief Construct a REECB controller from model data
       *
       * @param[in] bus Terminal bus the controller measures.
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(ReecbIdx::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::~Reecb()
      {
      }

      /**
       * @brief Terminal-bus voltage, real component
       *
       * @return Reference to the bus variable.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Reecb<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /**
       * @brief Terminal-bus voltage, imaginary component
       *
       * @return Reference to the bus variable.
       */
      template <typename scalar_type, typename index_type>
      scalar_type& Reecb<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      /**
       * @brief Resolve the parameter-derived constants and selector masks
       *
       * Raises each controller lag to the well-posedness floor, sizes the
       * component power base, and turns the four mode flags into complementary
       * multiplicative masks. The masks let the residual select control paths
       * without parameter-dependent control flow, which keeps its structure
       * fixed for sparse automatic differentiation.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::setDerivedParameters()
      {
        // The lags are raised to the floor in place, so a negative value is
        // rejected here while the value as read is still available. verify()
        // reports the count.
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Reecb: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(Trv_, "Trv");
        check_non_negative(Tp_, "Tp");
        check_non_negative(Tiq_, "Tiq");
        check_non_negative(Tpord_, "Tpord");

        if (Trv_ < TIME_CONSTANT_MINIMUM || Tp_ < TIME_CONSTANT_MINIMUM
            || Tiq_ < TIME_CONSTANT_MINIMUM || Tpord_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Reecb: Trv, Tp, Tiq, and Tpord below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the controller lags well posed\n";
        }

        Trv_   = std::max(Trv_, TIME_CONSTANT_MINIMUM);
        Tp_    = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tiq_   = std::max(Tiq_, TIME_CONSTANT_MINIMUM);
        Tpord_ = std::max(Tpord_, TIME_CONSTANT_MINIMUM);

        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);

        pf_on_ = PfFlag_ ? ONE<RealT> : ZERO<RealT>;
        v_on_  = VFlag_ ? ONE<RealT> : ZERO<RealT>;
        q_on_  = QFlag_ ? ONE<RealT> : ZERO<RealT>;

        pf_off_ = ONE<RealT> - pf_on_;
        v_off_  = ONE<RealT> - v_on_;
        q_off_  = ONE<RealT> - q_on_;

        p_priority_ = Pqflag_ ? ONE<RealT> : ZERO<RealT>;
        q_priority_ = ONE<RealT> - p_priority_;
      }

      /**
       * @brief Evaluate log(1 - exp(-x)) without cancellation
       *
       * Both terms approach one for a small argument, so the direct form loses
       * precision exactly where the limiter inversions below need it. The
       * hyperbolic identity is used under log 2 and log1p above it.
       *
       * @param[in] x Strictly positive argument.
       * @return The logarithm, always negative.
       */
      template <typename scalar_type, typename index_type>
      typename Reecb<scalar_type, index_type>::RealT
      Reecb<scalar_type, index_type>::logOneMinusExp(RealT x) const
      {
        static constexpr RealT LOG_TWO = static_cast<RealT>(0.6931471805599453);

        if (x < LOG_TWO)
        {
          return LOG_TWO - HALF<RealT> * x
                 + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      /**
       * @brief Recover the input that a smooth clamp maps to a requested output
       *
       * Initialization uses the same smooth CommonMath clamp as the residual,
       * so a steady state must be seeded with the limiter *input* rather than
       * its output. The smooth clamp is asymptotic at both limits, so a
       * requested output within the initialization tolerance of a limit is
       * represented by a finite offset past it instead of by the true infinite
       * preimage.
       *
       * @tparam LowerT Type of the lower limiter bound.
       * @tparam UpperT Type of the upper limiter bound.
       *
       * @param[in] requested_output Output the limiter must reproduce.
       * @param[in] lower_limit Lower limiter bound.
       * @param[in] upper_limit Upper limiter bound.
       * @param[out] limiter_input Input producing the requested output.
       * @return false when the request lies outside the limits, or when the
       *         limits coincide at a different value; the output is then unset.
       *
       * @note The limit types intentionally may differ from the scalar type so
       * that constant Real limits and algebraic-variable limits both work.
       */
      template <typename scalar_type, typename index_type>
      template <typename LowerT, typename UpperT>
      bool Reecb<scalar_type, index_type>::solveLimiterInput(
          ScalarT  requested_output,
          LowerT   lower_limit,
          UpperT   upper_limit,
          ScalarT& limiter_input) const
      {
        const RealT output_value = static_cast<RealT>(requested_output);
        const RealT lower_value  = static_cast<RealT>(lower_limit);
        const RealT upper_value  = static_cast<RealT>(upper_limit);

        if (lower_value > upper_value
            || output_value < lower_value - INITIALIZATION_TOLERANCE
            || output_value > upper_value + INITIALIZATION_TOLERANCE)
        {
          return false;
        }

        const RealT width = upper_value - lower_value;
        if (width <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(lower_value);
          return std::abs(output_value - lower_value) <= INITIALIZATION_TOLERANCE;
        }

        const RealT distance_from_lower = output_value - lower_value;
        const RealT distance_from_upper = upper_value - output_value;
        if (distance_from_lower <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(lower_value - INITIALIZATION_LIMIT_OFFSET);
          return true;
        }
        if (distance_from_upper <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(upper_value + INITIALIZATION_LIMIT_OFFSET);
          return true;
        }

        const RealT scaled_lower_distance = Math::MU<RealT> * distance_from_lower;
        const RealT scaled_upper_distance = Math::MU<RealT> * distance_from_upper;
        const RealT correction            = (scaled_lower_distance
                                  + logOneMinusExp(scaled_lower_distance)
                                  - logOneMinusExp(scaled_upper_distance))
                                 / Math::MU<RealT>;
        limiter_input = static_cast<ScalarT>(lower_value + correction);
        return true;
      }

      /**
       * @brief Choose a PI input whose anti-windup derivative is stationary
       *
       * An anti-windup integrator is at rest either because its rate is zero
       * or because the rate pushes into a limit that blocks it. A nonzero rate
       * is therefore parked just past the limit it drives toward.
       *
       * The zero-rate branch deliberately returns the nominal input without
       * clamping it: an inactive PI history may legitimately sit outside its
       * own output limits, and callers that need a representable output ask
       * solveLimiterInput() for one explicitly.
       *
       * @tparam LowerT Type of the lower limiter bound.
       * @tparam UpperT Type of the upper limiter bound.
       *
       * @param[in] nominal_input Input to keep when the rate is already zero.
       * @param[in] rate Anti-windup integrator rate.
       * @param[in] lower_limit Lower limiter bound.
       * @param[in] upper_limit Upper limiter bound.
       * @return A stationary integrator input.
       *
       * @note The limit types intentionally may differ from the scalar type so
       * that constant Real limits and algebraic-variable limits both work.
       */
      template <typename scalar_type, typename index_type>
      template <typename LowerT, typename UpperT>
      scalar_type Reecb<scalar_type, index_type>::steadyAntiWindupInput(
          ScalarT nominal_input,
          ScalarT rate,
          LowerT  lower_limit,
          UpperT  upper_limit) const
      {
        const RealT rate_value = static_cast<RealT>(rate);
        if (std::abs(rate_value) <= INITIALIZATION_TOLERANCE)
        {
          return nominal_input;
        }
        if (rate_value > ZERO<RealT>)
        {
          return upper_limit + static_cast<ScalarT>(INITIALIZATION_LIMIT_OFFSET);
        }
        return lower_limit - static_cast<ScalarT>(INITIALIZATION_LIMIT_OFFSET);
      }

      /**
       * @brief Convert a system-base power or current to REECB component base
       *
       * @param[in] value Quantity on the system base.
       * @return The same quantity on the component base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Reecb<scalar_type, index_type>::toComponentBase(
          scalar_type value) const
      {
        return value * va_system_base_ / va_converter_base_;
      }

      /**
       * @brief Convert a component-base power or current to the system base
       *
       * @param[in] value Quantity on the component base.
       * @return The same quantity on the system base.
       */
      template <typename scalar_type, typename index_type>
      scalar_type Reecb<scalar_type, index_type>::toSystemBase(
          scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      /**
       * @brief Read the parameters out of the model data
       *
       * Only the component power base is required; every other parameter keeps
       * the default documented in the model README when omitted. A missing
       * required key, a non-numeric value, or a switch outside {0, 1} is
       * counted and reported by verify() rather than throwing. Integer JSON
       * values are accepted for real parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;
        Vref0_given_           = false;

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
            Log::error() << "Reecb: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_switch = [&](auto key, bool& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* bool_value = std::get_if<bool>(&value))
          {
            target = *bool_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value);
                   index_value && (*index_value == 0 || *index_value == 1))
          {
            target = (*index_value == 1);
          }
          else
          {
            Log::error() << "Reecb: parameter '" << name
                         << "' must be bool or integer 0/1\n";
            ++parameter_error_count_;
          }
        };

        if (!data.parameters.contains(Params::mva))
        {
          Log::error() << "Reecb: missing required parameter 'mva'\n";
          ++parameter_error_count_;
        }
        load_real(Params::mva, mva_base_, "mva");
        load_switch(Params::PfFlag, PfFlag_, "PfFlag");
        load_switch(Params::VFlag, VFlag_, "VFlag");
        load_switch(Params::QFlag, QFlag_, "QFlag");
        load_switch(Params::Pqflag, Pqflag_, "Pqflag");
        load_real(Params::Trv, Trv_, "Trv");
        load_real(Params::Tp, Tp_, "Tp");
        if (data.parameters.contains(Params::Vref0))
        {
          load_real(Params::Vref0, Vref0_, "Vref0");
          Vref0_given_ = true;
        }
        load_real(Params::Vdip, Vdip_, "Vdip");
        load_real(Params::Vup, Vup_, "Vup");
        load_real(Params::dbd1, dbd1_, "dbd1");
        load_real(Params::dbd2, dbd2_, "dbd2");
        load_real(Params::kqv, kqv_, "kqv");
        load_real(Params::Iql1, Iql1_, "Iql1");
        load_real(Params::Iqh1, Iqh1_, "Iqh1");
        load_real(Params::Qmax, Qmax_, "Qmax");
        load_real(Params::Qmin, Qmin_, "Qmin");
        load_real(Params::Kqp, Kqp_, "Kqp");
        load_real(Params::Kqi, Kqi_, "Kqi");
        load_real(Params::Vmax, Vmax_, "Vmax");
        load_real(Params::Vmin, Vmin_, "Vmin");
        load_real(Params::Kvp, Kvp_, "Kvp");
        load_real(Params::Kvi, Kvi_, "Kvi");
        load_real(Params::Tiq, Tiq_, "Tiq");
        load_real(Params::Tpord, Tpord_, "Tpord");
        load_real(Params::dPmax, dPmax_, "dPmax");
        load_real(Params::dPmin, dPmin_, "dPmin");
        load_real(Params::Pmax, Pmax_, "Pmax");
        load_real(Params::Pmin, Pmin_, "Pmin");
        load_real(Params::Imax, Imax_, "Imax");
        setDerivedParameters();
      }

      /**
       * @brief Access the monitor
       *
       * @return Monitor for this model, or nullptr when the model was
       *         constructed without data.
       */
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Reecb<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Bind the monitorable variables to their internal states
       *
       * The two current commands are published on the system base and the two
       * filtered measurements on the component base, as documented in the
       * model README.
       */
      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeMonitor()
      {
        using I        = ReecbIdx;
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::iqcmd, [this]
                      { return y_.getData()[I::IQCMD]; });
        monitor_->set(Variable::ipcmd, [this]
                      { return y_.getData()[I::IPCMD]; });
        monitor_->set(Variable::vmeas, [this]
                      { return y_.getData()[I::VMEAS]; });
        monitor_->set(Variable::pmeas, [this]
                      { return y_.getData()[I::PMEAS]; });
      }

      /**
       * @brief Set the component ID
       *
       * @param[in] component_id Identifier assigned by the system model.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the model vectors and wire the command outputs
       *
       * Sizes the state, residual, bus-interface, and signal-interface
       * buffers, seeds the identity index maps, and points each assigned
       * command node at the internal state it publishes. Those nodes alias
       * REECB storage from here on, which is how initialize() reads the seeds
       * an upstream model wrote. Repeated calls reuse the allocated vectors.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::allocate()
      {
        using I = ReecbIdx;
        using E = ReecbExt;

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        auto signal_size = E::MAXIMUM;
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
              &y[I::IQCMD],
              &(this->getVariableIndex(static_cast<IdxT>(I::IQCMD))));
        }

        if (signals_.template isAssigned<ReecbInternalVariables::IPCMD>())
        {
          signals_.template getSignalNode<ReecbInternalVariables::IPCMD>()->set(
              &y[I::IPCMD],
              &(this->getVariableIndex(static_cast<IdxT>(I::IPCMD))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the REECB configuration
       *
       * Checks parameter-loading errors, static parameter relationships,
       * terminal-bus association, and attached external signals. Seeded
       * command feasibility is operating-point dependent and is checked by
       * initialize().
       *
       * @return int Number of configuration errors; zero when valid.
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

        if (bus_ == nullptr)
        {
          Log::error() << "Reecb: bus pointer is null\n";
          ret += 1;
        }

        check(mva_base_ > ZERO<RealT>, "mva must be positive");
        check(Vdip_ < Vup_, "Vdip must be less than Vup");
        check(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required");
        check(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
        check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
        check(Imax_ >= ZERO<RealT>, "Imax must be non-negative");

        // An attached port must resolve to writable signal storage. The
        // enumerator is a template argument, so each port names itself once.
        auto check_attached_signal =
            [&]<ReecbExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Reecb: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<ReecbExternalVariables::PE>("pe");
        check_attached_signal.template operator()<ReecbExternalVariables::QGEN>("qgen");
        check_attached_signal.template operator()<ReecbExternalVariables::QEXT>("qext");
        check_attached_signal.template operator()<ReecbExternalVariables::PFAREF>("pfaref");
        check_attached_signal.template operator()<ReecbExternalVariables::PREF>("pref");

        return ret;
      }

      /**
       * @brief Initialize REECB from seeded current-command ports
       *
       * Reads the assigned system-base `ipcmd` and `iqcmd` nodes, resolves a
       * component-base steady state that preserves those seeds, and initializes
       * attached feedback/reference signals. All operating-point checks are
       * completed before model or signal storage is modified.
       *
       * @pre allocate() has completed.
       * @pre verify() has reported no configuration errors.
       * @pre The terminal bus and assigned command nodes have been initialized.
       *
       * @return int 0 on success; nonzero when the commands are outside the
       *             current circle, the selected control path cannot represent
       *             them, or an initial reference is undefined.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initialize()
      {
        using I = ReecbIdx;

        auto* y = y_.getData();

        // Assigned command nodes alias these entries after allocate(). Their
        // system-base seeds remain untouched throughout initialization.
        const ScalarT ipcmd0_system = y[I::IPCMD];
        const ScalarT iqcmd0_system = y[I::IQCMD];
        const ScalarT ipcmd0        = toComponentBase(ipcmd0_system);
        const ScalarT iqcmd0        = toComponentBase(iqcmd0_system);
        const RealT   ipcmd0_value  = static_cast<RealT>(ipcmd0);
        const RealT   iqcmd0_value  = static_cast<RealT>(iqcmd0);

        const ScalarT vr          = Vr();
        const ScalarT vi          = Vi();
        const ScalarT vt0         = std::sqrt(vr * vr + vi * vi);
        const ScalarT vmeas0      = vt0;
        const ScalarT vmeas_safe0 = Math::max(vmeas0, VMEAS_MINIMUM);
        const ScalarT pmeas0      = ipcmd0 * vmeas_safe0;
        const ScalarT qgen0       = iqcmd0 * vmeas_safe0;
        const RealT   vref0       = Vref0_given_ ? Vref0_ : static_cast<RealT>(vt0);

        if (!std::isfinite(ipcmd0_value) || !std::isfinite(iqcmd0_value) || !std::isfinite(static_cast<RealT>(vt0)))
        {
          Log::error() << "Reecb: initial bus voltage and current commands must be finite\n";
          return 1;
        }
        if (ipcmd0_value < ZERO<RealT>)
        {
          Log::error() << "Reecb: initial active-current command must be non-negative\n";
          return 1;
        }

        const RealT current_squared0       = ipcmd0_value * ipcmd0_value + iqcmd0_value * iqcmd0_value;
        const RealT current_limit_squared0 = Imax_ * Imax_;
        if (current_squared0 > current_limit_squared0 + INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Reecb: initial current commands exceed the Imax circle\n";
          return 1;
        }

        const RealT pmeas0_value = static_cast<RealT>(pmeas0);
        if (pmeas0_value < Pmin_ - INITIALIZATION_TOLERANCE || pmeas0_value > Pmax_ + INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Reecb: initial active power is outside Pmin/Pmax\n";
          return 1;
        }

        const RealT iqcirc_squared0 = current_limit_squared0 - p_priority_ * ipcmd0_value * ipcmd0_value;
        const RealT ipcirc_squared0 = current_limit_squared0 - q_priority_ * iqcmd0_value * iqcmd0_value;
        if (iqcirc_squared0 < -INITIALIZATION_TOLERANCE || ipcirc_squared0 < -INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Reecb: initial current commands violate the selected priority circle\n";
          return 1;
        }

        const ScalarT iqcirc0 = static_cast<ScalarT>(std::sqrt(std::max(iqcirc_squared0, ZERO<RealT>)));
        const ScalarT ipcirc0 = static_cast<ScalarT>(std::sqrt(std::max(ipcirc_squared0, ZERO<RealT>)));
        const ScalarT iqmax0  = q_priority_ * static_cast<ScalarT>(Imax_) + p_priority_ * iqcirc0;
        const ScalarT ipmax0  = p_priority_ * static_cast<ScalarT>(Imax_) + q_priority_ * ipcirc0;

        const ScalarT sdip0  = Math::inside(vt0, Vdip_, Vup_);
        const ScalarT verr0  = Math::deadband2(static_cast<ScalarT>(vref0) - vmeas0, dbd1_, dbd2_);
        const ScalarT iqv0   = Math::clamp(kqv_ * verr0, Iql1_, Iqh1_);
        const ScalarT iqinj0 = (ONE<RealT> - sdip0) * iqv0;

        ScalarT iqraw0{};
        if (!solveLimiterInput(iqcmd0, -iqmax0, iqmax0, iqraw0))
        {
          Log::error() << "Reecb: initial reactive-current command is outside the available current limit\n";
          return 1;
        }

        ScalarT ip_limiter_input0{};
        if (!solveLimiterInput(ipcmd0, ZERO<RealT>, ipmax0, ip_limiter_input0))
        {
          Log::error() << "Reecb: initial active-current command is outside the available current limit\n";
          return 1;
        }
        const ScalarT pord0 = ip_limiter_input0 * vmeas_safe0;

        ScalarT fpord0{};
        if (!solveLimiterInput(static_cast<ScalarT>(ZERO<RealT>), dPmin_, dPmax_, fpord0))
        {
          Log::error() << "Reecb: zero initial active-power ramp is outside dPmin/dPmax\n";
          return 1;
        }
        const ScalarT rpord0 = Math::clamp(fpord0, dPmin_, dPmax_);
        const ScalarT pref0  = pord0 + Tpord_ * fpord0;

        const ScalarT iq_control0 = iqraw0 - iqinj0;

        struct ReactiveSeed
        {
          ScalarT qv{};
          ScalarT qref{};
          ScalarT eq{};
          ScalarT vpiq{};
          ScalarT epiv{};
          ScalarT xpiq{};
          ScalarT iqbase{};
          ScalarT xpiv{};
        } reactive;

        if (!QFlag_)
        {
          reactive.qv   = iq_control0;
          reactive.qref = reactive.qv * vmeas_safe0;
        }
        else if (!VFlag_)
        {
          reactive.qref = vmeas0;
          reactive.qv   = reactive.qref / vmeas_safe0;
        }
        else
        {
          if (!solveLimiterInput(qgen0, Qmin_, Qmax_, reactive.qref))
          {
            Log::error() << "Reecb: initial reactive power is outside Qmin/Qmax\n";
            return 1;
          }

          reactive.qv = reactive.qref / vmeas_safe0;
        }

        reactive.eq = Math::clamp(reactive.qref, Qmin_, Qmax_) - qgen0;

        ScalarT vpiq_input0{};
        if (QFlag_ && VFlag_)
        {
          if (!solveLimiterInput(vmeas0, Vmin_, Vmax_, vpiq_input0))
          {
            Log::error() << "Reecb: initial voltage is outside Vmin/Vmax\n";
            return 1;
          }
        }
        else
        {
          const ScalarT vpiq_nominal0 = v_on_ * vmeas0 + v_off_ * reactive.qref;
          vpiq_input0                 = steadyAntiWindupInput(vpiq_nominal0, Kqi_ * reactive.eq, Vmin_, Vmax_);
        }

        reactive.vpiq = Math::clamp(vpiq_input0, Vmin_, Vmax_);
        reactive.epiv = v_on_ * reactive.vpiq + v_off_ * reactive.qref - vmeas0;
        reactive.xpiq = vpiq_input0 - Kqp_ * reactive.eq;

        ScalarT iqbase_input0{};
        if (QFlag_)
        {
          if (!solveLimiterInput(iq_control0, -iqmax0, iqmax0, iqbase_input0))
          {
            Log::error() << "Reecb: initial reactive-current command is outside the voltage-controller current limit\n";
            return 1;
          }
        }
        else
        {
          iqbase_input0 = steadyAntiWindupInput(static_cast<ScalarT>(ZERO<RealT>), Kvi_ * reactive.epiv, -iqmax0, iqmax0);
        }

        reactive.iqbase = Math::clamp(iqbase_input0, -iqmax0, iqmax0);
        reactive.xpiv   = iqbase_input0 - Kvp_ * reactive.epiv;

        ScalarT pfaref0 = static_cast<ScalarT>(ZERO<RealT>);
        if (PfFlag_)
        {
          if (std::abs(pmeas0_value) <= INITIALIZATION_TOLERANCE)
          {
            if (std::abs(static_cast<RealT>(reactive.qref)) > INITIALIZATION_TOLERANCE)
            {
              Log::error() << "Reecb: power-factor control cannot represent nonzero Qref at zero active power\n";
              return 1;
            }
          }
          else
          {
            pfaref0 = static_cast<ScalarT>(std::atan(static_cast<RealT>(reactive.qref / pmeas0)));
          }
        }

        const ScalarT pe0_system   = toSystemBase(pmeas0);
        const ScalarT qgen0_system = toSystemBase(qgen0);
        const ScalarT qext0_system = toSystemBase(reactive.qref);
        const ScalarT pref0_system = toSystemBase(pref0);

        y[I::VMEAS]     = vmeas0;
        y[I::PMEAS]     = pmeas0;
        y[I::XPIQ]      = reactive.xpiq;
        y[I::XPIV]      = reactive.xpiv;
        y[I::QV]        = reactive.qv;
        y[I::PORD]      = pord0;
        y[I::VT]        = vt0;
        y[I::VMEASSAFE] = vmeas_safe0;
        y[I::SDIP]      = sdip0;
        y[I::VERR]      = verr0;
        y[I::IQV]       = iqv0;
        y[I::QREF]      = reactive.qref;
        y[I::EQ]        = reactive.eq;
        y[I::VPIQ]      = reactive.vpiq;
        y[I::EPIV]      = reactive.epiv;
        y[I::FPORD]     = fpord0;
        y[I::RPORD]     = rpord0;
        y[I::IQCIRC]    = iqcirc0;
        y[I::IPCIRC]    = ipcirc0;
        y[I::IQMAX]     = iqmax0;
        y[I::IPMAX]     = ipmax0;
        y[I::IQBASE]    = reactive.iqbase;
        y[I::IQRAW]     = iqraw0;

        if (!Vref0_given_)
        {
          Vref0_ = vref0;
        }

        pe_set_     = pe0_system;
        qgen_set_   = qgen0_system;
        qext_set_   = qext0_system;
        pfaref_set_ = pfaref0;
        pref_set_   = pref0_system;

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::PE>(pe_set_);
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          signals_.template writeExternalVariable<ReecbExternalVariables::QGEN>(qgen_set_);
        }
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
       * The two measurement filters, the two PI states, the reactive-current
       * lag, and the active-power order carry derivatives; every other
       * internal variable is algebraic.
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::tagDifferentiable()
      {
        using I = ReecbIdx;

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[I::VMEAS] = true;
        tag_[I::PMEAS] = true;
        tag_[I::XPIQ]  = true;
        tag_[I::XPIV]  = true;
        tag_[I::QV]    = true;
        tag_[I::PORD]  = true;
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * All REECB variables are per-unit voltages, powers, or currents of the
       * same order, so they share the relative tolerance as their absolute
       * floor.
       *
       * @param[in] rel_tol Solver relative tolerance.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Internal residual
       *
       * Evaluates the six controller states and the nineteen algebraic rows
       * documented in the model README. The body is kept free of branches and
       * loops so that sparse automatic differentiation resolves a fixed
       * structure; the four mode selections enter as the multiplicative masks
       * set by setDerivedParameters().
       *
       * @param[in] y Internal variables.
       * @param[in] yp Internal variable derivatives.
       * @param[in] wb Terminal-bus voltage components.
       * @param[in] ws External signal values on system base.
       * @param[out] f Internal residuals.
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Reecb<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        using I = ReecbIdx;
        using E = ReecbExt;

        const ScalarT vmeas        = y[I::VMEAS];
        const ScalarT pmeas        = y[I::PMEAS];
        const ScalarT xpiq         = y[I::XPIQ];
        const ScalarT xpiv         = y[I::XPIV];
        const ScalarT qv           = y[I::QV];
        const ScalarT pord         = y[I::PORD];
        const ScalarT vt           = y[I::VT];
        const ScalarT vmeas_safe   = y[I::VMEASSAFE];
        const ScalarT sdip         = y[I::SDIP];
        const ScalarT verr         = y[I::VERR];
        const ScalarT iqv          = y[I::IQV];
        const ScalarT qref         = y[I::QREF];
        const ScalarT eq           = y[I::EQ];
        const ScalarT vpiq         = y[I::VPIQ];
        const ScalarT epiv         = y[I::EPIV];
        const ScalarT fpord        = y[I::FPORD];
        const ScalarT rpord        = y[I::RPORD];
        const ScalarT iqcirc       = y[I::IQCIRC];
        const ScalarT ipcirc       = y[I::IPCIRC];
        const ScalarT iqmax        = y[I::IQMAX];
        const ScalarT ipmax        = y[I::IPMAX];
        const ScalarT iqbase       = y[I::IQBASE];
        const ScalarT iqraw        = y[I::IQRAW];
        const ScalarT iqcmd_system = y[I::IQCMD];
        const ScalarT ipcmd_system = y[I::IPCMD];

        const ScalarT vmeas_dot = yp[I::VMEAS];
        const ScalarT pmeas_dot = yp[I::PMEAS];
        const ScalarT xpiq_dot  = yp[I::XPIQ];
        const ScalarT xpiv_dot  = yp[I::XPIV];
        const ScalarT qv_dot    = yp[I::QV];
        const ScalarT pord_dot  = yp[I::PORD];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT pe     = toComponentBase(ws[E::PE]);
        const ScalarT qgen   = toComponentBase(ws[E::QGEN]);
        const ScalarT qext   = toComponentBase(ws[E::QEXT]);
        const ScalarT pfaref = ws[E::PFAREF];
        const ScalarT pref   = toComponentBase(ws[E::PREF]);
        const ScalarT iqcmd  = toComponentBase(iqcmd_system);
        const ScalarT ipcmd  = toComponentBase(ipcmd_system);

        f[I::VMEAS]     = -vmeas_dot + (vt - vmeas) / Trv_;
        f[I::PMEAS]     = -pmeas_dot + (pe - pmeas) / Tp_;
        f[I::XPIQ]      = -xpiq_dot + sdip * Math::antiwindup(Kqp_ * eq + xpiq, Kqi_ * eq, Vmin_, Vmax_);
        f[I::XPIV]      = -xpiv_dot + sdip * Math::antiwindup(Kvp_ * epiv + xpiv, Kvi_ * epiv, -iqmax, iqmax);
        f[I::QV]        = -qv_dot + sdip * (qref / vmeas_safe - qv) / Tiq_;
        f[I::PORD]      = -pord_dot + sdip * Math::antiwindup(pord, rpord, Pmin_, Pmax_);
        f[I::VT]        = -vt * vt + vr * vr + vi * vi;
        f[I::VMEASSAFE] = -vmeas_safe + Math::max(vmeas, VMEAS_MINIMUM);
        f[I::SDIP]      = -sdip + Math::inside(vt, Vdip_, Vup_);
        f[I::VERR]      = -verr + Math::deadband2(Vref0_ - vmeas, dbd1_, dbd2_);
        f[I::IQV]       = -iqv + Math::clamp(kqv_ * verr, Iql1_, Iqh1_);
        f[I::QREF]      = -qref + pf_on_ * pmeas * std::tan(pfaref) + pf_off_ * qext;
        f[I::EQ]        = -eq + Math::clamp(qref, Qmin_, Qmax_) - qgen;
        f[I::VPIQ]      = -vpiq + Math::clamp(Kqp_ * eq + xpiq, Vmin_, Vmax_);
        f[I::EPIV]      = -epiv + v_on_ * vpiq + v_off_ * qref - vmeas;
        f[I::FPORD]     = -fpord + (pref - pord) / Tpord_;
        f[I::RPORD]     = -rpord + Math::clamp(fpord, dPmin_, dPmax_);
        f[I::IQCIRC]    = -iqcirc * iqcirc + Imax_ * Imax_ - p_priority_ * ipcmd * ipcmd;
        f[I::IPCIRC]    = -ipcirc * ipcirc + Imax_ * Imax_ - q_priority_ * iqcmd * iqcmd;
        f[I::IQMAX]     = -iqmax + q_priority_ * Imax_ + p_priority_ * iqcirc;
        f[I::IPMAX]     = -ipmax + p_priority_ * Imax_ + q_priority_ * ipcirc;
        f[I::IQBASE]    = -iqbase + Math::clamp(Kvp_ * epiv + xpiv, -iqmax, iqmax);
        f[I::IQRAW]     = -iqraw + q_on_ * iqbase + q_off_ * qv + (ONE<RealT> - sdip) * iqv;
        f[I::IQCMD]     = -iqcmd_system + toSystemBase(Math::clamp(iqraw, -iqmax, iqmax));
        f[I::IPCMD]     = -ipcmd_system + toSystemBase(Math::clamp(pord / vmeas_safe, ZERO<RealT>, ipmax));

        return 0;
      }

      /**
       * @brief Residuals of system equations
       *
       * Refreshes the bus and signal interface buffers and evaluates the
       * internal residual. REECB injects no current, so there is no bus
       * residual. An unattached input port falls back to the value latched by
       * initialize().
       *
       * @return int 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::evaluateResidual()
      {
        using E = ReecbExt;

        ws_[E::PE]     = pe_set_;
        ws_[E::QGEN]   = qgen_set_;
        ws_[E::QEXT]   = qext_set_;
        ws_[E::PFAREF] = pfaref_set_;
        ws_[E::PREF]   = pref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<ReecbExternalVariables::PE>())
        {
          ws_[E::PE]         = signals_.template readExternalVariable<ReecbExternalVariables::PE>();
          ws_indices_[E::PE] = signals_.template readExternalVariableIndex<ReecbExternalVariables::PE>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QGEN>())
        {
          ws_[E::QGEN]         = signals_.template readExternalVariable<ReecbExternalVariables::QGEN>();
          ws_indices_[E::QGEN] = signals_.template readExternalVariableIndex<ReecbExternalVariables::QGEN>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::QEXT>())
        {
          ws_[E::QEXT]         = signals_.template readExternalVariable<ReecbExternalVariables::QEXT>();
          ws_indices_[E::QEXT] = signals_.template readExternalVariableIndex<ReecbExternalVariables::QEXT>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PFAREF>())
        {
          ws_[E::PFAREF]         = signals_.template readExternalVariable<ReecbExternalVariables::PFAREF>();
          ws_indices_[E::PFAREF] = signals_.template readExternalVariableIndex<ReecbExternalVariables::PFAREF>();
        }
        if (signals_.template isAttached<ReecbExternalVariables::PREF>())
        {
          ws_[E::PREF]         = signals_.template readExternalVariable<ReecbExternalVariables::PREF>();
          ws_indices_[E::PREF] = signals_.template readExternalVariableIndex<ReecbExternalVariables::PREF>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();
        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
