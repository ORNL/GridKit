/**
 * @file RepcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REPCA phasor-dynamics plant-control model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REPCA/RepcaData.hpp>
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
       * @brief Construct REPCA with its documented parameter defaults.
       * @param[in] bus Regulated bus measured by the controller.
       */
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RepcaIdx::MAXIMUM);
        setDerivedParameters();
      }

      /**
       * @brief Construct REPCA from parameters and monitor selections.
       * @param[in] bus Regulated bus measured by the controller.
       * @param[in] data Model parameters and monitor selections.
       */
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RepcaIdx::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::~Repca()
      {
      }

      /// Regulated-bus voltage, real component.
      template <typename scalar_type, typename index_type>
      scalar_type& Repca<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /// Regulated-bus voltage, imaginary component.
      template <typename scalar_type, typename index_type>
      scalar_type& Repca<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      /**
       * @brief Load provided parameters and retain documented defaults for omissions.
       *
       * Numeric parameters accept real or integer values. Mode switches accept
       * bool and numeric 0/1 values; invalid types are counted for verify().
       */
      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
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
            Log::error() << "Repca: parameter '" << name << "' must be numeric\n";
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
          else if (const auto* real_value = std::get_if<RealT>(&value);
                   real_value && (*real_value == RealT{0} || *real_value == RealT{1}))
          {
            target = (*real_value == RealT{1});
          }
          else
          {
            Log::error() << "Repca: parameter '" << name << "' must be bool or 0/1\n";
            ++parameter_error_count_;
          }
        };

        load_real(Params::mva, mva_base_, "mva");
        load_switch(Params::VcompFlag, VcompFlag_, "VcompFlag");
        load_switch(Params::RefFlag, RefFlag_, "RefFlag");
        load_switch(Params::Freqflag, Freqflag_, "Freqflag");
        load_real(Params::Tfltr, Tfltr_, "Tfltr");
        load_real(Params::Vfrz, Vfrz_, "Vfrz");
        load_real(Params::Rc, Rc_, "Rc");
        load_real(Params::Xc, Xc_, "Xc");
        load_real(Params::Kc, Kc_, "Kc");
        load_real(Params::dbdlow, dbdlow_, "dbdlow");
        load_real(Params::dbdupper, dbdupper_, "dbdupper");
        load_real(Params::emax, emax_, "emax");
        load_real(Params::emin, emin_, "emin");
        load_real(Params::Kp, Kp_, "Kp");
        load_real(Params::Ki, Ki_, "Ki");
        load_real(Params::Qmax, Qmax_, "Qmax");
        load_real(Params::Qmin, Qmin_, "Qmin");
        load_real(Params::Tft, Tft_, "Tft");
        load_real(Params::Tfv, Tfv_, "Tfv");
        load_real(Params::Tp, Tp_, "Tp");
        load_real(Params::fdbd1, fdbd1_, "fdbd1");
        load_real(Params::fdbd2, fdbd2_, "fdbd2");
        load_real(Params::Ddn, Ddn_, "Ddn");
        load_real(Params::Dup, Dup_, "Dup");
        load_real(Params::femax, femax_, "femax");
        load_real(Params::femin, femin_, "femin");
        load_real(Params::Kpg, Kpg_, "Kpg");
        load_real(Params::Kig, Kig_, "Kig");
        load_real(Params::Pmax, Pmax_, "Pmax");
        load_real(Params::Pmin, Pmin_, "Pmin");
        load_real(Params::Tlag, Tlag_, "Tlag");
        setDerivedParameters();
      }

      /**
       * @brief Resolve time floors, the component power base, and mode masks.
       */
      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::setDerivedParameters()
      {
        if (Tfltr_ < TIME_CONSTANT_MINIMUM || Tp_ < TIME_CONSTANT_MINIMUM
            || Tlag_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Repca: Tfltr, Tp, and Tlag below " << TIME_CONSTANT_MINIMUM
                         << " s are raised to that floor to keep the controller lags well posed\n";
        }

        Tfltr_ = std::max(Tfltr_, TIME_CONSTANT_MINIMUM);
        Tp_    = std::max(Tp_, TIME_CONSTANT_MINIMUM);
        Tlag_  = std::max(Tlag_, TIME_CONSTANT_MINIMUM);

        va_converter_base_ = mva_base_ * static_cast<RealT>(1.0e6);

        vcomp_on_  = VcompFlag_ ? ONE<RealT> : ZERO<RealT>;
        vcomp_off_ = ONE<RealT> - vcomp_on_;
        ref_on_    = RefFlag_ ? ONE<RealT> : ZERO<RealT>;
        ref_off_   = ONE<RealT> - ref_on_;
        freq_on_   = Freqflag_ ? ONE<RealT> : ZERO<RealT>;
      }

      /// Evaluate log(1 - exp(-x)) accurately for a positive argument.
      template <typename scalar_type, typename index_type>
      typename Repca<scalar_type, index_type>::RealT
      Repca<scalar_type, index_type>::logOneMinusExp(RealT x) const
      {
        const RealT log_two = std::log(static_cast<RealT>(2.0));

        if (x < log_two)
        {
          return log_two - HALF<RealT> * x + std::log(std::sinh(HALF<RealT> * x));
        }
        return std::log1p(-std::exp(-x));
      }

      /**
       * @brief Recover an input for a requested smooth-clamp output.
       * @return false when the limits are invalid or the requested output lies outside them.
       */
      template <typename scalar_type, typename index_type>
      bool Repca<scalar_type, index_type>::solveLimiterInput(
          ScalarT  requested_output,
          RealT    lower_limit,
          RealT    upper_limit,
          ScalarT& limiter_input) const
      {
        const RealT output_value = static_cast<RealT>(requested_output);

        if (lower_limit > upper_limit || output_value < lower_limit || output_value > upper_limit)
        {
          return false;
        }

        const RealT width = upper_limit - lower_limit;
        if (width <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(lower_limit);
          return true;
        }

        const RealT distance_from_lower = output_value - lower_limit;
        const RealT distance_from_upper = upper_limit - output_value;
        if (distance_from_lower <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(lower_limit - INITIALIZATION_LIMIT_OFFSET);
          return true;
        }
        if (distance_from_upper <= INITIALIZATION_TOLERANCE)
        {
          limiter_input = static_cast<ScalarT>(upper_limit + INITIALIZATION_LIMIT_OFFSET);
          return true;
        }

        const RealT mu                    = Math::MU<RealT>;
        const RealT scaled_lower_distance = mu * distance_from_lower;
        const RealT scaled_upper_distance = mu * distance_from_upper;
        const RealT log_lower             = logOneMinusExp(scaled_lower_distance);
        const RealT log_upper             = logOneMinusExp(scaled_upper_distance);
        const RealT correction            = (scaled_lower_distance + log_lower - log_upper) / mu;

        limiter_input = static_cast<ScalarT>(lower_limit + correction);
        return true;
      }

      /// Convert a system-base power quantity to REPCA component base.
      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_converter_base_;
      }

      /// Convert a component-base power quantity to system base.
      template <typename scalar_type, typename index_type>
      scalar_type Repca<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      /// Access the configured monitor, or nullptr for the default constructor.
      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Repca<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /// Bind monitor selections to REPCA state entries in public monitor order.
      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeMonitor()
      {
        using I        = RepcaIdx;
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::qext, [this]
                      { return y_.getData()[I::QEXT]; });
        monitor_->set(Variable::pext, [this]
                      { return y_.getData()[I::PEXT]; });
        monitor_->set(Variable::vmeas, [this]
                      { return y_.getData()[I::VMEAS]; });
        monitor_->set(Variable::qmeas, [this]
                      { return y_.getData()[I::QMEAS]; });
        monitor_->set(Variable::pmeas, [this]
                      { return y_.getData()[I::PMEAS]; });
      }

      /// Set the component identifier used by monitor labels.
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate model buffers and bind assigned command-output nodes.
       * @return 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::allocate()
      {
        using I = RepcaIdx;
        using E = RepcaExt;

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        const auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});

        const auto signal_size = E::MAXIMUM;
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        auto* y = y_.getData();

        if (signals_.template isAssigned<RepcaInternalVariables::QEXT>())
        {
          signals_.template getSignalNode<RepcaInternalVariables::QEXT>()->set(
              &y[I::QEXT],
              &(this->getVariableIndex(static_cast<IdxT>(I::QEXT))));
        }

        if (signals_.template isAssigned<RepcaInternalVariables::PEXT>())
        {
          signals_.template getSignalNode<RepcaInternalVariables::PEXT>()->set(
              &y[I::PEXT],
              &(this->getVariableIndex(static_cast<IdxT>(I::PEXT))));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate parameters, the regulated bus, and signal connections.
       * @return The number of configuration errors.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Repca: " << message << '\n';
            ret += 1;
          }
        };

        if (bus_ == nullptr)
        {
          Log::error() << "Repca: bus pointer is null\n";
          ret += 1;
        }

        check(mva_base_ > ZERO<RealT>, "mva must be positive");
        check(Tfv_ > ZERO<RealT>, "Tfv must be positive");
        check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_,
              "dbdlow <= 0 <= dbdupper is required");
        check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_,
              "emin <= 0 <= emax is required");
        check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
        check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_,
              "fdbd1 <= 0 <= fdbd2 is required");
        check(femin_ <= ZERO<RealT> && ZERO<RealT> <= femax_,
              "femin <= 0 <= femax is required");
        check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");

        auto check_required_signal = [&]<RepcaExternalVariables variable>(const char* name)
        {
          if (!signals_.template isAttached<variable>())
          {
            Log::error() << "Repca: " << name << " signal is required\n";
            ret += 1;
          }
          else if (!signals_.template isLinked<variable>())
          {
            Log::error() << "Repca: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_required_signal.template operator()<RepcaExternalVariables::IBRANCHR>("ibranchr");
        check_required_signal.template operator()<RepcaExternalVariables::IBRANCHI>("ibranchi");
        check_required_signal.template operator()<RepcaExternalVariables::PBRANCH>("pbranch");
        check_required_signal.template operator()<RepcaExternalVariables::QBRANCH>("qbranch");
        check_required_signal.template operator()<RepcaExternalVariables::FREQ>("freq");

        auto check_optional_signal = [&]<RepcaExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "Repca: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_optional_signal.template operator()<RepcaExternalVariables::FREQREF>("freqref");
        check_optional_signal.template operator()<RepcaExternalVariables::VREF>("vref");
        check_optional_signal.template operator()<RepcaExternalVariables::QREF>("qref");
        check_optional_signal.template operator()<RepcaExternalVariables::PPLANTREF>("pplantref");

        return ret;
      }

      /**
       * @brief Initialize from seeded qext and, when frequency control is enabled, pext.
       *
       * All feasibility checks precede state, derivative, latch, and signal
       * writes, so failure leaves model and signal storage unchanged.
       *
       * @pre allocate() has completed and verify() returned zero.
       * @return 0 on success; nonzero for an inadmissible operating point.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize()
      {
        using I = RepcaIdx;

        auto* y = y_.getData();

        const ScalarT qext0_system = y[I::QEXT];
        const ScalarT pext0_system = y[I::PEXT];
        const ScalarT qext0        = toComponentBase(qext0_system);
        const ScalarT pext0        = toComponentBase(pext0_system);

        const ScalarT vr = Vr();
        const ScalarT vi = Vi();
        const ScalarT ibranchr =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHR>();
        const ScalarT ibranchi =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHI>();
        const ScalarT pbranch_system =
            signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
        const ScalarT qbranch_system =
            signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
        const ScalarT freq = signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();

        auto is_finite = [](ScalarT value)
        {
          return std::isfinite(static_cast<RealT>(value));
        };
        if (!is_finite(vr) || !is_finite(vi) || !is_finite(ibranchr)
            || !is_finite(ibranchi) || !is_finite(pbranch_system)
            || !is_finite(qbranch_system)
            || !is_finite(freq) || !is_finite(qext0)
            || (Freqflag_ && !is_finite(pext0)))
        {
          Log::error() << "Repca: initial bus, signal, and command values must be finite\n";
          return 1;
        }

        const ScalarT pbranch = toComponentBase(pbranch_system);
        const ScalarT qbranch = toComponentBase(qbranch_system);

        const ScalarT vldc_r = vr - Rc_ * ibranchr + Xc_ * ibranchi;
        const ScalarT vldc_i = vi - Rc_ * ibranchi - Xc_ * ibranchr;

        const ScalarT v0      = std::sqrt(vr * vr + vi * vi);
        const ScalarT vldc0   = std::sqrt(vldc_r * vldc_r + vldc_i * vldc_i);
        const ScalarT vdroop0 = v0 + Kc_ * qbranch;
        const ScalarT vctrl0  = vcomp_on_ * vldc0 + vcomp_off_ * vdroop0;
        const ScalarT vmeas0  = vctrl0;
        const ScalarT qmeas0  = qbranch;
        const ScalarT pmeas0  = pbranch;
        const ScalarT sfrz0   = Math::above(v0, Vfrz_);

        const ScalarT zero    = static_cast<ScalarT>(ZERO<RealT>);
        const ScalarT erq0    = zero;
        const ScalarT erqdb0  = Math::deadband2(erq0, dbdlow_, dbdupper_);
        const ScalarT erqlim0 = Math::clamp(erqdb0, emin_, emax_);
        const ScalarT qpi0    = qext0;
        const ScalarT xqlag0  = qpi0;

        ScalarT qpi_input0{};
        if (!solveLimiterInput(qpi0, Qmin_, Qmax_, qpi_input0))
        {
          Log::error() << "Repca: initial reactive-power command is outside Qmin/Qmax\n";
          return 1;
        }
        const ScalarT xqpi0      = qpi_input0 - Kp_ * erqlim0;
        const ScalarT q_aw_rate0 = Math::antiwindup(qpi0, Ki_ * erqlim0, Qmin_, Qmax_);
        if (!is_finite(q_aw_rate0)
            || std::abs(static_cast<RealT>(q_aw_rate0)) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Repca: reactive-power PI antiwindup rate is nonzero at initialization\n";
          return 1;
        }

        const ScalarT ef0    = Math::deadband2(zero, fdbd1_, fdbd2_);
        const ScalarT pfreq0 = Ddn_ * Math::ramp(ef0) - Dup_ * Math::ramp(-ef0);
        const ScalarT ep0    = zero;
        const ScalarT eplim0 = Math::clamp(ep0, femin_, femax_);
        const ScalarT pref0  = Freqflag_ ? pext0 : Math::clamp(pmeas0, Pmin_, Pmax_);
        const ScalarT ppi0   = pref0;
        ScalarT       ppi_input0{};
        if (Freqflag_)
        {
          if (!solveLimiterInput(ppi0, Pmin_, Pmax_, ppi_input0))
          {
            Log::error() << "Repca: initial active-power command is outside Pmin/Pmax\n";
            return 1;
          }
        }
        else
        {
          ppi_input0 = pmeas0;
        }
        const ScalarT xppi0      = ppi_input0 - Kpg_ * eplim0;
        const ScalarT p_aw_rate0 = Math::antiwindup(ppi0, Kig_ * eplim0, Pmin_, Pmax_);
        if (!is_finite(p_aw_rate0)
            || std::abs(static_cast<RealT>(p_aw_rate0)) > INITIALIZATION_TOLERANCE)
        {
          Log::error() << "Repca: active-power PI antiwindup rate is nonzero at initialization\n";
          return 1;
        }

        const ScalarT pext_output0      = Freqflag_ ? pext0_system : zero;
        const ScalarT freqref0          = freq;
        const ScalarT vref0             = vmeas0;
        const ScalarT qref0_system      = qbranch_system;
        const ScalarT pplantref0_system = toSystemBase(pmeas0 - pfreq0);

        y[I::VMEAS]  = vmeas0;
        y[I::QMEAS]  = qmeas0;
        y[I::XQPI]   = xqpi0;
        y[I::XQLAG]  = xqlag0;
        y[I::PMEAS]  = pmeas0;
        y[I::XPPI]   = xppi0;
        y[I::PREF]   = pref0;
        y[I::V]      = v0;
        y[I::VLDC]   = vldc0;
        y[I::VDROOP] = vdroop0;
        y[I::VCTRL]  = vctrl0;
        y[I::SFRZ]   = sfrz0;
        y[I::ERQ]    = erq0;
        y[I::ERQDB]  = erqdb0;
        y[I::ERQLIM] = erqlim0;
        y[I::QPI]    = qpi0;
        y[I::QEXT]   = qext0_system;
        y[I::EF]     = ef0;
        y[I::EP]     = ep0;
        y[I::EPLIM]  = eplim0;
        y[I::PPI]    = ppi0;
        y[I::PEXT]   = pext_output0;

        freqref_set_   = freqref0;
        vref_set_      = vref0;
        qref_set_      = qref0_system;
        pplantref_set_ = pplantref0_system;

        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::FREQREF>(
              freqref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::VREF>(vref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::QREF>(qref_set_);
        }
        if (signals_.template isAttached<RepcaExternalVariables::PPLANTREF>())
        {
          signals_.template writeExternalVariable<RepcaExternalVariables::PPLANTREF>(
              pplantref_set_);
        }

        y_.setDataUpdated();
        yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));
        return 0;
      }

      /**
       * @brief Tag the seven controller states as differential variables.
       * @return 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::tagDifferentiable()
      {
        using I = RepcaIdx;

        std::fill(tag_.begin(), tag_.end(), false);
        tag_[I::VMEAS] = true;
        tag_[I::QMEAS] = true;
        tag_[I::XQPI]  = true;
        tag_[I::XQLAG] = true;
        tag_[I::PMEAS] = true;
        tag_[I::XPPI]  = true;
        tag_[I::PREF]  = true;
        return 0;
      }

      /**
       * @brief Set the common absolute-tolerance floor for all REPCA variables.
       * @return 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Evaluate the seven differential and fifteen algebraic rows.
       *
       * Mode selections enter through fixed masks so the residual remains
       * branch-free for sparse automatic differentiation.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Repca<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        using I = RepcaIdx;
        using E = RepcaExt;

        const ScalarT vmeas  = y[I::VMEAS];
        const ScalarT qmeas  = y[I::QMEAS];
        const ScalarT xqpi   = y[I::XQPI];
        const ScalarT xqlag  = y[I::XQLAG];
        const ScalarT pmeas  = y[I::PMEAS];
        const ScalarT xppi   = y[I::XPPI];
        const ScalarT pref   = y[I::PREF];
        const ScalarT v      = y[I::V];
        const ScalarT vldc   = y[I::VLDC];
        const ScalarT vdroop = y[I::VDROOP];
        const ScalarT vctrl  = y[I::VCTRL];
        const ScalarT sfrz   = y[I::SFRZ];
        const ScalarT erq    = y[I::ERQ];
        const ScalarT erqdb  = y[I::ERQDB];
        const ScalarT erqlim = y[I::ERQLIM];
        const ScalarT qpi    = y[I::QPI];
        const ScalarT qext   = toComponentBase(y[I::QEXT]);
        const ScalarT ef     = y[I::EF];
        const ScalarT ep     = y[I::EP];
        const ScalarT eplim  = y[I::EPLIM];
        const ScalarT ppi    = y[I::PPI];
        const ScalarT pext   = toComponentBase(y[I::PEXT]);

        const ScalarT vmeas_dot = yp[I::VMEAS];
        const ScalarT qmeas_dot = yp[I::QMEAS];
        const ScalarT xqpi_dot  = yp[I::XQPI];
        const ScalarT xqlag_dot = yp[I::XQLAG];
        const ScalarT pmeas_dot = yp[I::PMEAS];
        const ScalarT xppi_dot  = yp[I::XPPI];
        const ScalarT pref_dot  = yp[I::PREF];

        const ScalarT vr = wb[0];
        const ScalarT vi = wb[1];

        const ScalarT ibranchr  = ws[E::IBRANCHR];
        const ScalarT ibranchi  = ws[E::IBRANCHI];
        const ScalarT pbranch   = toComponentBase(ws[E::PBRANCH]);
        const ScalarT qbranch   = toComponentBase(ws[E::QBRANCH]);
        const ScalarT freq      = ws[E::FREQ];
        const ScalarT freqref   = ws[E::FREQREF];
        const ScalarT vref      = ws[E::VREF];
        const ScalarT qref      = toComponentBase(ws[E::QREF]);
        const ScalarT pplantref = toComponentBase(ws[E::PPLANTREF]);

        const ScalarT vldc_r = vr - Rc_ * ibranchr + Xc_ * ibranchi;
        const ScalarT vldc_i = vi - Rc_ * ibranchi - Xc_ * ibranchr;
        const ScalarT pfreq  = Ddn_ * Math::ramp(ef) - Dup_ * Math::ramp(-ef);

        f[I::VMEAS] = -vmeas_dot + (vctrl - vmeas) / Tfltr_;
        f[I::QMEAS] = -qmeas_dot + (qbranch - qmeas) / Tfltr_;
        f[I::XQPI]  = -xqpi_dot + sfrz * Math::antiwindup(qpi, Ki_ * erqlim, Qmin_, Qmax_);
        f[I::XQLAG] = -xqlag_dot + (qpi - xqlag) / Tfv_;
        f[I::PMEAS] = -pmeas_dot + (pbranch - pmeas) / Tp_;
        f[I::XPPI]  = -xppi_dot + Math::antiwindup(ppi, Kig_ * eplim, Pmin_, Pmax_);
        f[I::PREF]  = -pref_dot + (ppi - pref) / Tlag_;

        f[I::V]      = -v * v + vr * vr + vi * vi;
        f[I::VLDC]   = -vldc * vldc + vldc_r * vldc_r + vldc_i * vldc_i;
        f[I::VDROOP] = -vdroop + v + Kc_ * qbranch;
        f[I::VCTRL]  = -vctrl + vcomp_on_ * vldc + vcomp_off_ * vdroop;
        f[I::SFRZ]   = -sfrz + Math::above(v, Vfrz_);
        f[I::ERQ]    = -erq + ref_on_ * (vref - vmeas) + ref_off_ * (qref - qmeas);
        f[I::ERQDB]  = -erqdb + Math::deadband2(erq, dbdlow_, dbdupper_);
        f[I::ERQLIM] = -erqlim + Math::clamp(erqdb, emin_, emax_);
        f[I::QPI]    = -qpi + Math::clamp(Kp_ * erqlim + xqpi, Qmin_, Qmax_);
        f[I::QEXT]   = -Tfv_ * (qext - xqlag) + Tft_ * (qpi - xqlag);

        f[I::EF]    = -ef + Math::deadband2(freqref - freq, fdbd1_, fdbd2_);
        f[I::EP]    = -ep + pplantref - pmeas + pfreq;
        f[I::EPLIM] = -eplim + Math::clamp(ep, femin_, femax_);
        f[I::PPI]   = -ppi + Math::clamp(Kpg_ * eplim + xppi, Pmin_, Pmax_);
        f[I::PEXT]  = -pext + freq_on_ * pref;

        return 0;
      }

      /**
       * @brief Refresh bus and signal buffers and evaluate the REPCA residual.
       *
       * Required measurements are read directly; unattached optional references
       * use the values latched by initialize(). REPCA has no bus residual.
       *
       * @return 0 on success.
       */
      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::evaluateResidual()
      {
        using E = RepcaExt;

        ws_[E::FREQREF]   = freqref_set_;
        ws_[E::VREF]      = vref_set_;
        ws_[E::QREF]      = qref_set_;
        ws_[E::PPLANTREF] = pplantref_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        ws_[E::IBRANCHR] =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHR>();
        ws_indices_[E::IBRANCHR] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::IBRANCHR>();
        ws_[E::IBRANCHI] =
            signals_.template readExternalVariable<RepcaExternalVariables::IBRANCHI>();
        ws_indices_[E::IBRANCHI] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::IBRANCHI>();
        ws_[E::PBRANCH] = signals_.template readExternalVariable<RepcaExternalVariables::PBRANCH>();
        ws_indices_[E::PBRANCH] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::PBRANCH>();
        ws_[E::QBRANCH] = signals_.template readExternalVariable<RepcaExternalVariables::QBRANCH>();
        ws_indices_[E::QBRANCH] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::QBRANCH>();
        ws_[E::FREQ] = signals_.template readExternalVariable<RepcaExternalVariables::FREQ>();
        ws_indices_[E::FREQ] =
            signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQ>();

        if (signals_.template isAttached<RepcaExternalVariables::FREQREF>())
        {
          ws_[E::FREQREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::FREQREF>();
          ws_indices_[E::FREQREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::FREQREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::VREF>())
        {
          ws_[E::VREF] = signals_.template readExternalVariable<RepcaExternalVariables::VREF>();
          ws_indices_[E::VREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::VREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::QREF>())
        {
          ws_[E::QREF] = signals_.template readExternalVariable<RepcaExternalVariables::QREF>();
          ws_indices_[E::QREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::QREF>();
        }
        if (signals_.template isAttached<RepcaExternalVariables::PPLANTREF>())
        {
          ws_[E::PPLANTREF] =
              signals_.template readExternalVariable<RepcaExternalVariables::PPLANTREF>();
          ws_indices_[E::PPLANTREF] =
              signals_.template readExternalVariableIndex<RepcaExternalVariables::PPLANTREF>();
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
