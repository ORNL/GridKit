#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Load/PERC1/Perc1.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    Perc1<scalar_type, index_type>::Perc1(BusT* bus)
      : bus_(bus), parameter_error_count_(2)
    {
      size_ = static_cast<IdxT>(Utilities::enum_size<Perc1InternalVariables>());
    }

    template <typename scalar_type, typename index_type>
    Perc1<scalar_type, index_type>::Perc1(BusT* bus, const ModelDataT& data)
      : bus_(bus),
        monitor_(std::make_unique<MonitorT>(data))
    {
      size_ = static_cast<IdxT>(Utilities::enum_size<Perc1InternalVariables>());
      initializeParameters(data);
      initializeMonitor();
    }

    template <typename scalar_type, typename index_type>
    Perc1<scalar_type, index_type>::~Perc1() = default;

    /// Load numeric parameters without throwing; verify() reports loading errors.
    template <typename scalar_type, typename index_type>
    void Perc1<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = typename ModelDataT::Parameters;

      auto load_real = [&](Params key, RealT& target, const char* name, bool required = false)
      {
        if (!data.parameters.contains(key))
        {
          if (required)
          {
            Log::error() << "Perc1: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
          }
          return;
        }

        const auto& value = data.parameters.at(key);
        RealT       parsed{};
        if (const auto* real_value = std::get_if<RealT>(&value))
        {
          parsed = *real_value;
        }
        else if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          parsed = static_cast<RealT>(*index_value);
        }
        else
        {
          Log::error() << "Perc1: parameter '" << name << "' must be numeric\n";
          ++parameter_error_count_;
          return;
        }
        if (!std::isfinite(parsed))
        {
          Log::error() << "Perc1: parameter '" << name << "' must be finite\n";
          ++parameter_error_count_;
          return;
        }
        target = parsed;
      };

      load_real(Params::Pnom, Pnom_, "Pnom", true);
      load_real(Params::Qnom, Qnom_, "Qnom", true);
      load_real(Params::Lfm, Lfm_, "Lfm");
      load_real(Params::QPratio, QPratio_, "QPratio");
      load_real(Params::Dbfl, Dbfl_, "Dbfl");
      load_real(Params::Dbfh, Dbfh_, "Dbfh");
      load_real(Params::Kdroop, Kdroop_, "Kdroop");
      load_real(Params::Kvp, Kvp_, "Kvp");
      load_real(Params::Tvp, Tvp_, "Tvp");
      load_real(Params::Kvq, Kvq_, "Kvq");
      load_real(Params::Tvq, Tvq_, "Tvq");
      load_real(Params::Tap, Tap_, "Tap");
      load_real(Params::Tbp, Tbp_, "Tbp");
      load_real(Params::Taq, Taq_, "Taq");
      load_real(Params::Tbq, Tbq_, "Tbq");
      load_real(Params::nP, nP_, "nP");
      load_real(Params::nQ, nQ_, "nQ");
      load_real(Params::Ipmax, Ipmax_, "Ipmax");
      load_real(Params::Ipmin, Ipmin_, "Ipmin");
      load_real(Params::Iqmax, Iqmax_, "Iqmax");
      load_real(Params::Iqmin, Iqmin_, "Iqmin");
      load_real(Params::Fcease, Fcease_, "Fcease");
      load_real(Params::Vcease, Vcease_, "Vcease");
      load_real(Params::Tcease, Tcease_, "Tcease");
      load_real(Params::Tdelay, Tdelay_, "Tdelay");
      load_real(Params::Vrecon, Vrecon_, "Vrecon");
      load_real(Params::Trecon, Trecon_, "Trecon");
      load_real(Params::Tramp, Tramp_, "Tramp");
      load_real(Params::Frecon, Frecon_, "Frecon");
      load_real(Params::Tt, Tt_, "Tt");
      load_real(Params::Tv, Tv_, "Tv");
      load_real(Params::Tf, Tf_, "Tf");
    }

    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      const auto size = static_cast<size_t>(size_);
      tag_.assign(size, false);
      variable_indices_.resize(size);
      residual_indices_.resize(size);
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }
      wb_.resize(2);
      h_.resize(2);
      const auto signal_size = Utilities::enum_size<Perc1ExternalVariables>();
      ws_.resize(static_cast<IdxT>(signal_size));
      ws_.setToZero();
      ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);
      allocated_ = true;
      return 0;
    }

    /// Validate the raw parameters before flooring denominator time constants.
    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::verify() const
    {
      int  ret   = parameter_error_count_;
      auto check = [&](bool condition, const char* message)
      {
        if (!condition)
        {
          Log::error() << "Perc1: " << message << '\n';
          ++ret;
        }
      };

      check(bus_ != nullptr, "terminal bus is required");
      check(std::isfinite(va_system_base_) && va_system_base_ > ZERO<RealT>,
            "system power base must be finite and positive");
      check(Pnom_ >= ZERO<RealT>, "Pnom must be non-negative");
      check(Dbfl_ <= ZERO<RealT> && Dbfh_ >= ZERO<RealT>, "Dbfl <= 0 <= Dbfh is required");
      check(Tvp_ >= ZERO<RealT> && Tvq_ >= ZERO<RealT>
                && Tap_ >= ZERO<RealT> && Tbp_ >= ZERO<RealT>
                && Taq_ >= ZERO<RealT> && Tbq_ >= ZERO<RealT>
                && Tt_ >= ZERO<RealT> && Tv_ >= ZERO<RealT> && Tf_ >= ZERO<RealT>,
            "electrical time constants must be non-negative");
      check(Ipmin_ <= Ipmax_, "Ipmin must not exceed Ipmax");
      check(Iqmin_ <= Iqmax_, "Iqmin must not exceed Iqmax");
      check(ZERO<RealT> <= Fcease_ && Fcease_ <= ONE<RealT>, "Fcease must lie in [0, 1]");
      check(Frecon_ >= ZERO<RealT>, "Frecon must be non-negative");
      check(ZERO<RealT> <= Vcease_ && Vcease_ <= Vrecon_, "0 <= Vcease <= Vrecon is required");
      check(Tcease_ >= ZERO<RealT> && Tdelay_ >= ZERO<RealT>
                && Trecon_ >= ZERO<RealT> && Tramp_ >= ZERO<RealT>,
            "cessation and reconnection durations must be non-negative");
      const auto freq = ports_.in.template port<Perc1SignalInputs::freq>();
      check(!freq.connected() || freq.linked(), "freq signal attached with no linked source");
      return ret;
    }

    template <typename scalar_type, typename index_type>
    void Perc1<scalar_type, index_type>::setDerivedParameters()
    {
      if (Tbp_ < TIME_CONSTANT_MINIMUM || Tbq_ < TIME_CONSTANT_MINIMUM
          || Tvp_ < TIME_CONSTANT_MINIMUM || Tvq_ < TIME_CONSTANT_MINIMUM
          || Tv_ < TIME_CONSTANT_MINIMUM || Tf_ < TIME_CONSTANT_MINIMUM || Tt_ < TIME_CONSTANT_MINIMUM)
      {
        static std::atomic<bool> warned{false};
        if (!warned.exchange(true))
        {
          Log::warning() << "Perc1: electrical denominator time constants below "
                         << TIME_CONSTANT_MINIMUM << " s are raised to preserve Hessenberg form\n";
        }
      }
      Tbp_ = std::max(Tbp_, TIME_CONSTANT_MINIMUM);
      Tbq_ = std::max(Tbq_, TIME_CONSTANT_MINIMUM);
      Tvp_ = std::max(Tvp_, TIME_CONSTANT_MINIMUM);
      Tvq_ = std::max(Tvq_, TIME_CONSTANT_MINIMUM);
      Tv_  = std::max(Tv_, TIME_CONSTANT_MINIMUM);
      Tf_  = std::max(Tf_, TIME_CONSTANT_MINIMUM);
      Tt_  = std::max(Tt_, TIME_CONSTANT_MINIMUM);
      if (Lfm_ < static_cast<RealT>(0.001))
      {
        Lfm_ = static_cast<RealT>(0.80);
      }
      auto timer = [](RealT duration, RealT& enabled, RealT& inverse)
      {
        enabled = duration > ZERO<RealT> ? ONE<RealT> : ZERO<RealT>;
        inverse = duration > ZERO<RealT> ? ONE<RealT> / duration : ZERO<RealT>;
      };
      timer(Tcease_, timed_cease_, inv_Tcease_);
      timer(Tdelay_, timed_delay_, inv_Tdelay_);
      timer(Trecon_, timed_recon_, inv_Trecon_);
      timer(Tramp_, timed_ramp_, inv_Tramp_);
    }

    template <typename scalar_type, typename index_type>
    typename Perc1<scalar_type, index_type>::RealT
    Perc1<scalar_type, index_type>::inverseClamp(RealT value, RealT lower, RealT upper)
    {
      // Equality is intentional: a collapsed interval has one feasible target.
      if (lower == upper && value == lower)
      {
        return value;
      }
      if (!(lower < value && value < upper))
      {
        return std::numeric_limits<RealT>::quiet_NaN();
      }
      return value + (std::log(-std::expm1(-Math::MU<RealT> * (value - lower))) - std::log(-std::expm1(-Math::MU<RealT> * (upper - value)))) / Math::MU<RealT>;
    }

    /**
     * @brief Preserve the initial bus voltage, frequency, and consumed P/Q.
     *
     * The constant-impedance fallback is selected once per initialization.
     * History starts empty with derivatives taken from the smooth equations.
     * @pre allocate() and terminal-bus initialization have completed.
     */
    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::initialize()
    {
      if (verify() != 0)
      {
        return 1;
      }

      // Initialization uses real values so references retain no AD dependencies.
      const RealT vr   = static_cast<RealT>(bus_->Vr());
      const RealT vi   = static_cast<RealT>(bus_->Vi());
      RealT       freq = ONE<RealT>;
      if (auto port = ports_.in.template port<Perc1SignalInputs::freq>())
      {
        freq = static_cast<RealT>(port.readSignal());
      }
      const RealT vnom2     = vr * vr + vi * vi;
      const bool  zero_load = Pnom_ == ZERO<RealT> && Qnom_ == ZERO<RealT>;
      if (!std::isfinite(vr) || !std::isfinite(vi) || !std::isfinite(freq)
          || !std::isfinite(vnom2) || (!zero_load && !(vnom2 > ZERO<RealT>) ))
      {
        Log::error() << "Perc1: initial inputs must be finite and a nonzero load requires nonzero voltage\n";
        return 1;
      }

      const bool               static_load = zero_load || Pnom_ == ZERO<RealT> || std::sqrt(vnom2) < Vcease_;
      constexpr auto           count       = Utilities::enum_size<Perc1InternalVariables>();
      std::array<RealT, count> initial{};
      // Local aliases follow the same variable order as the residual.
      const auto               XPLAG  = static_cast<size_t>(Perc1InternalVariables::XPLAG);
      const auto               XQLAG  = static_cast<size_t>(Perc1InternalVariables::XQLAG);
      const auto               XPWO   = static_cast<size_t>(Perc1InternalVariables::XPWO);
      const auto               XQWO   = static_cast<size_t>(Perc1InternalVariables::XQWO);
      const auto               DFMEAS = static_cast<size_t>(Perc1InternalVariables::DFMEAS);
      const auto               VMEAS  = static_cast<size_t>(Perc1InternalVariables::VMEAS);
      const auto               IP     = static_cast<size_t>(Perc1InternalVariables::IP);
      const auto               IQ     = static_cast<size_t>(Perc1InternalVariables::IQ);
      const auto               XD     = static_cast<size_t>(Perc1InternalVariables::XD);
      const auto               TA     = static_cast<size_t>(Perc1InternalVariables::TA);
      const auto               V      = static_cast<size_t>(Perc1InternalVariables::V);
      const auto               VSAFE  = static_cast<size_t>(Perc1InternalVariables::VSAFE);
      const auto               UP     = static_cast<size_t>(Perc1InternalVariables::UP);
      const auto               UQ     = static_cast<size_t>(Perc1InternalVariables::UQ);
      const auto               PCTRL  = static_cast<size_t>(Perc1InternalVariables::PCTRL);
      const auto               QCTRL  = static_cast<size_t>(Perc1InternalVariables::QCTRL);
      const auto               IPV    = static_cast<size_t>(Perc1InternalVariables::IPV);
      const auto               IQV    = static_cast<size_t>(Perc1InternalVariables::IQV);
      const auto               FRACON = static_cast<size_t>(Perc1InternalVariables::FRACON);
      const auto               IR     = static_cast<size_t>(Perc1InternalVariables::IR);
      const auto               II     = static_cast<size_t>(Perc1InternalVariables::II);

      initial[FRACON] = ONE<RealT>;
      RealT g         = ZERO<RealT>;
      RealT b         = ZERO<RealT>;
      RealT pref      = ZERO<RealT>;
      RealT qref      = ZERO<RealT>;
      RealT iqextra   = ZERO<RealT>;
      RealT vnom_safe = ONE<RealT>;
      if (!zero_load)
      {
        g           = Pnom_ / vnom2;
        b           = Qnom_ / vnom2;
        initial[IR] = -g * vr - b * vi;
        initial[II] = -g * vi + b * vr;
      }
      setDerivedParameters();
      const RealT derived[] = {Tap_ / Tbp_, Taq_ / Tbq_, Kvp_ / Tvp_, Kvq_ / Tvq_, inv_Tcease_, inv_Tdelay_, inv_Trecon_, inv_Tramp_};
      for (const RealT value : derived)
      {
        if (!std::isfinite(value))
        {
          Log::error() << "Perc1: derived gains and timer reciprocals must be finite\n";
          return 1;
        }
      }
      const RealT v = std::sqrt(vnom2 + VOLTAGE_EPSILON * VOLTAGE_EPSILON);
      vnom_safe     = Math::max(v, VOLTAGE_MINIMUM);
      initial[V] = initial[VMEAS] = initial[XPWO] = initial[XQWO] = v;
      initial[VSAFE]                                              = vnom_safe;
      if (!static_load)
      {
        this->setComponentBase(va_system_base_ * Pnom_ / Lfm_);
        const RealT kbase         = va_system_base_ / va_component_base_;
        const RealT inverse_kbase = va_component_base_ / va_system_base_;
        if (!std::isfinite(va_component_base_) || !(va_component_base_ > ZERO<RealT>)
            || !std::isfinite(kbase) || !(kbase > ZERO<RealT>)
            || !std::isfinite(inverse_kbase) || !(inverse_kbase > ZERO<RealT>) )
        {
          Log::error() << "Perc1: component power base and conversion ratios must be finite and positive\n";
          return 1;
        }

        initial[DFMEAS] = Math::above(v, static_cast<RealT>(0.7)) * (freq - ONE<RealT>);

        const RealT s0    = Math::sigmoid(-HALF<RealT>);
        const RealT ramp0 = Tramp_ > ZERO<RealT> ? Math::clamp(ZERO<RealT>, ZERO<RealT>, ONE<RealT>) : ONE<RealT>;
        initial[FRACON]   = ONE<RealT> - Fcease_
                          + Fcease_ * ((ONE<RealT> - s0) + s0 * (ONE<RealT> - s0) * Frecon_ * ramp0);
        initial[IP]  = kbase * Pnom_ * v / vnom2;
        initial[IQ]  = kbase * Qnom_ * v / vnom2;
        initial[IPV] = initial[IP] / initial[FRACON];
        initial[IQV] = QPratio_ * initial[IPV];
        iqextra      = initial[IQ] / initial[FRACON] - initial[IQV];

        const RealT up = vnom_safe * inverseClamp(initial[IPV], Ipmin_, Ipmax_);
        const RealT uq = vnom_safe * inverseClamp(initial[IQV], Iqmin_, Iqmax_);
        if (!std::isfinite(up) || !std::isfinite(uq))
        {
          Log::error() << "Perc1: initial controlled currents must lie inside noncollapsed limits or equal a collapsed bound\n";
          return 1;
        }
        initial[UP] = initial[XPLAG] = initial[PCTRL] = up;
        initial[UQ] = initial[XQLAG] = initial[QCTRL] = uq;
        pref                                          = up - Kdroop_ * Math::deadband2(initial[DFMEAS], Dbfl_, Dbfh_);
        qref                                          = uq;
      }

      if (!std::all_of(initial.begin(), initial.end(), [](RealT value)
                       { return std::isfinite(value); })
          || !std::isfinite(g) || !std::isfinite(b) || !std::isfinite(pref)
          || !std::isfinite(qref) || !std::isfinite(iqextra) || !std::isfinite(vnom_safe))
      {
        Log::error() << "Perc1: initialization produced a nonfinite value\n";
        return 1;
      }

      dynamic_on_ = static_load ? ZERO<RealT> : ONE<RealT>;
      static_on_  = ONE<RealT> - dynamic_on_;
      cease_on_   = !static_load && Fcease_ > ZERO<RealT> ? ONE<RealT> : ZERO<RealT>;
      G_          = static_on_ * g;
      B_          = static_on_ * b;
      if (static_load)
      {
        this->setComponentBase(va_system_base_);
      }
      Pref_      = pref;
      Qref_      = qref;
      Iqextra_   = iqextra;
      Vnom_safe_ = vnom_safe;

      auto* y  = y_.getData();
      auto* yp = yp_.getData();
      for (size_t i = 0; i < count; ++i)
      {
        y[i]  = initial[i];
        yp[i] = ZERO<RealT>;
      }
      // Evaluate only the internal equations: initialization must not stamp the bus.
      const ScalarT              wb[] = {static_cast<ScalarT>(vr), static_cast<ScalarT>(vi)};
      const ScalarT              ws[] = {static_cast<ScalarT>(freq)};
      std::array<ScalarT, count> residual{};
      evaluateInternalResidual(y, yp, wb, ws, residual.data());
      for (size_t i = XD; i <= TA; ++i)
      {
        yp[i] = static_cast<RealT>(residual[i]);
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
      {
        this->initializeDependencyTrackingVariableNumbers();
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::tagDifferentiable()
    {
      std::fill(tag_.begin(), tag_.end(), false);
      std::fill_n(tag_.begin(), static_cast<size_t>(Perc1InternalVariables::V), true);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Evaluate the documented smooth PERC1 equations.
     *
     * Configuration masks and safe timer reciprocals keep the expression
     * branchless for sparse AD, including static loads and zero-duration timers.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int
    Perc1<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT* y, const ScalarT* yp, const ScalarT* wb, const ScalarT* ws, ScalarT* f)
    {
      const auto    XPLAG  = static_cast<size_t>(Perc1InternalVariables::XPLAG);
      const auto    XQLAG  = static_cast<size_t>(Perc1InternalVariables::XQLAG);
      const auto    XPWO   = static_cast<size_t>(Perc1InternalVariables::XPWO);
      const auto    XQWO   = static_cast<size_t>(Perc1InternalVariables::XQWO);
      const auto    DFMEAS = static_cast<size_t>(Perc1InternalVariables::DFMEAS);
      const auto    VMEAS  = static_cast<size_t>(Perc1InternalVariables::VMEAS);
      const auto    IP     = static_cast<size_t>(Perc1InternalVariables::IP);
      const auto    IQ     = static_cast<size_t>(Perc1InternalVariables::IQ);
      const auto    XD     = static_cast<size_t>(Perc1InternalVariables::XD);
      const auto    XC     = static_cast<size_t>(Perc1InternalVariables::XC);
      const auto    XR     = static_cast<size_t>(Perc1InternalVariables::XR);
      const auto    XH     = static_cast<size_t>(Perc1InternalVariables::XH);
      const auto    TC     = static_cast<size_t>(Perc1InternalVariables::TC);
      const auto    TD     = static_cast<size_t>(Perc1InternalVariables::TD);
      const auto    TR     = static_cast<size_t>(Perc1InternalVariables::TR);
      const auto    TA     = static_cast<size_t>(Perc1InternalVariables::TA);
      const auto    V      = static_cast<size_t>(Perc1InternalVariables::V);
      const auto    VSAFE  = static_cast<size_t>(Perc1InternalVariables::VSAFE);
      const auto    DP     = static_cast<size_t>(Perc1InternalVariables::DP);
      const auto    DQ     = static_cast<size_t>(Perc1InternalVariables::DQ);
      const auto    UP     = static_cast<size_t>(Perc1InternalVariables::UP);
      const auto    UQ     = static_cast<size_t>(Perc1InternalVariables::UQ);
      const auto    PCTRL  = static_cast<size_t>(Perc1InternalVariables::PCTRL);
      const auto    QCTRL  = static_cast<size_t>(Perc1InternalVariables::QCTRL);
      const auto    IPV    = static_cast<size_t>(Perc1InternalVariables::IPV);
      const auto    IQV    = static_cast<size_t>(Perc1InternalVariables::IQV);
      const auto    FRACON = static_cast<size_t>(Perc1InternalVariables::FRACON);
      const auto    IR     = static_cast<size_t>(Perc1InternalVariables::IR);
      const auto    II     = static_cast<size_t>(Perc1InternalVariables::II);
      const auto    FREQ   = static_cast<size_t>(Perc1ExternalVariables::FREQ);
      const ScalarT vr     = wb[0];
      const ScalarT vi     = wb[1];

      const ScalarT xplag  = y[XPLAG];
      const ScalarT xqlag  = y[XQLAG];
      const ScalarT xpwo   = y[XPWO];
      const ScalarT xqwo   = y[XQWO];
      const ScalarT dfmeas = y[DFMEAS];
      const ScalarT vmeas  = y[VMEAS];
      const ScalarT ip     = y[IP];
      const ScalarT iq     = y[IQ];
      const ScalarT xd     = y[XD];
      const ScalarT xc     = y[XC];
      const ScalarT xr     = y[XR];
      const ScalarT xh     = y[XH];
      const ScalarT tc     = y[TC];
      const ScalarT td     = y[TD];
      const ScalarT tr     = y[TR];
      const ScalarT ta     = y[TA];
      const ScalarT v      = y[V];
      const ScalarT vsafe  = y[VSAFE];
      const ScalarT dp     = y[DP];
      const ScalarT dq     = y[DQ];
      const ScalarT up     = y[UP];
      const ScalarT uq     = y[UQ];
      const ScalarT pctrl  = y[PCTRL];
      const ScalarT qctrl  = y[QCTRL];
      const ScalarT ipv    = y[IPV];
      const ScalarT iqv    = y[IQV];
      const ScalarT fracon = y[FRACON];
      const ScalarT ir     = y[IR];
      const ScalarT ii     = y[II];

      const ScalarT sd        = Math::sigmoid(xd - HALF<RealT>);
      const ScalarT sc        = Math::sigmoid(xc - HALF<RealT>);
      const ScalarT sr        = Math::sigmoid(xr - HALF<RealT>);
      const ScalarT sh        = Math::sigmoid(xh - HALF<RealT>);
      const ScalarT cease_run = (ONE<RealT> - sc)
                                * (ONE<RealT> - Math::sigmoid(xd - static_cast<RealT>(0.75)))
                                * Math::sigmoid(Vcease_ - vmeas);
      const ScalarT recon_run = sc * (ONE<RealT> - sd) * sh
                                * Math::sigmoid(vmeas - Vrecon_) * Math::sigmoid(ONE<RealT> - ta / TAU);
      const ScalarT cease_done = ONE<RealT> - timed_cease_ + timed_cease_ * Math::sigmoid(tc * inv_Tcease_ - ONE<RealT>);
      const ScalarT delay_done = ONE<RealT> - timed_delay_ + timed_delay_ * Math::sigmoid(td * inv_Tdelay_ - ONE<RealT>);
      const ScalarT recon_done = ONE<RealT> - timed_recon_ + timed_recon_ * Math::sigmoid(tr * inv_Trecon_ - ONE<RealT>);
      const ScalarT ramp       = ONE<RealT> - timed_ramp_ + timed_ramp_ * Math::clamp(ta * inv_Tramp_, ZERO<RealT>, ONE<RealT>);

      f[XPLAG]  = -yp[XPLAG] + dynamic_on_ * (up - xplag) / Tbp_;
      f[XQLAG]  = -yp[XQLAG] + dynamic_on_ * (uq - xqlag) / Tbq_;
      f[XPWO]   = -yp[XPWO] + dynamic_on_ * (v - xpwo) / Tvp_;
      f[XQWO]   = -yp[XQWO] + dynamic_on_ * (v - xqwo) / Tvq_;
      f[DFMEAS] = -yp[DFMEAS] + dynamic_on_ * (Math::above(v, static_cast<RealT>(0.7)) * (ws[FREQ] - ONE<RealT>) -dfmeas) / Tf_;
      f[VMEAS]  = -yp[VMEAS] + dynamic_on_ * (v - vmeas) / Tv_;
      f[IP]     = -yp[IP] + dynamic_on_ * (fracon * ipv - ip) / Tt_;
      f[IQ]     = -yp[IQ] + dynamic_on_ * (fracon * (iqv + Iqextra_) - iq) / Tt_;
      f[XD]     = -yp[XD] + cease_on_ * latch(xd, cease_run * cease_done, sc) / TAU;
      f[XC]     = -yp[XC] + cease_on_ * latch(xc, sd * delay_done, sr) / TAU;
      f[XR]     = -yp[XR] + cease_on_ * latch(xr, recon_run * recon_done, sd) / TAU;
      f[XH]     = -yp[XH] + cease_on_ * latch(xh, sc, static_cast<ScalarT>(ZERO<RealT>)) / TAU;
      f[TC]     = -yp[TC] + cease_on_ * (timed_cease_ * cease_run - (ONE<RealT> - timed_cease_ * cease_run) * tc / TAU);
      f[TD]     = -yp[TD] + cease_on_ * (timed_delay_ * sd - (ONE<RealT> - timed_delay_ * sd) * td / TAU);
      f[TR]     = -yp[TR] + cease_on_ * (timed_recon_ * recon_run - (ONE<RealT> - timed_recon_ * recon_run) * tr / TAU);
      f[TA]     = -yp[TA] + cease_on_ * (timed_ramp_ * sr * (ONE<RealT> - sc) * (ONE<RealT> - sd) * (ONE<RealT> - Math::sigmoid(ta * inv_Tramp_ - ONE<RealT>)) - (ONE<RealT> - timed_ramp_ + timed_ramp_ * sc) * ta / TAU);

      f[V]      = -v * v + vr * vr + vi * vi + VOLTAGE_EPSILON * VOLTAGE_EPSILON;
      f[VSAFE]  = -vsafe + Math::max(vmeas, VOLTAGE_MINIMUM);
      f[DP]     = -dp + dynamic_on_ * Kvp_ * (v - xpwo) / Tvp_;
      f[DQ]     = -dq + dynamic_on_ * Kvq_ * (v - xqwo) / Tvq_;
      f[UP]     = -up + Pref_ + dp + dynamic_on_ * Kdroop_ * Math::deadband2(dfmeas, Dbfl_, Dbfh_);
      f[UQ]     = -uq + Qref_ + dq;
      f[PCTRL]  = -pctrl + xplag + (Tap_ / Tbp_) * (up - xplag);
      f[QCTRL]  = -qctrl + xqlag + (Taq_ / Tbq_) * (uq - xqlag);
      f[IPV]    = -ipv + dynamic_on_ * Math::clamp(pctrl / vsafe * std::exp(nP_ * std::log(vsafe / Vnom_safe_)), Ipmin_, Ipmax_);
      f[IQV]    = -iqv + dynamic_on_ * Math::clamp(qctrl / vsafe * std::exp(nQ_ * std::log(vsafe / Vnom_safe_)), Iqmin_, Iqmax_);
      f[FRACON] = -fracon + ONE<RealT> - dynamic_on_ * Fcease_
                  + dynamic_on_ * Fcease_ * ((ONE<RealT> - sh) + sh * (ONE<RealT> - sc) * Frecon_ * ramp);
      f[IR] = -dynamic_on_ * (this->toComponentBase(v * ir) + vr * ip + vi * iq)
              - static_on_ * ir - G_ * vr - B_ * vi;
      f[II] = -dynamic_on_ * (this->toComponentBase(v * ii) + vi * ip - vr * iq)
              - static_on_ * ii - G_ * vi + B_ * vr;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int
    Perc1<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT* y, [[maybe_unused]] const ScalarT* yp, [[maybe_unused]] const ScalarT* wb, ScalarT* h)
    {
      h[0] = y[static_cast<size_t>(Perc1InternalVariables::IR)];
      h[1] = y[static_cast<size_t>(Perc1InternalVariables::II)];
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::evaluateResidual()
    {
      const auto FREQ   = static_cast<size_t>(Perc1ExternalVariables::FREQ);
      auto*      ws     = ws_.getData();
      ws[FREQ]          = ONE<RealT>;
      ws_indices_[FREQ] = INVALID_INDEX<IdxT>;
      if (auto port = ports_.in.template port<Perc1SignalInputs::freq>())
      {
        ws[FREQ]          = port.readSignal();
        ws_indices_[FREQ] = port.signalVariableIndex();
      }
      auto* wb = wb_.getData();
      wb[0]    = bus_->Vr();
      wb[1]    = bus_->Vi();
      evaluateInternalResidual(y_.getData(), yp_.getData(), wb, ws, f_.getData());
      evaluateBusResidual(y_.getData(), yp_.getData(), wb, h_.getData());
      bus_->Ir() += h_.getData()[0];
      bus_->Ii() += h_.getData()[1];
      if (bus_->size() > 0)
      {
        bus_->getResidual().setDataUpdated();
      }
      f_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Perc1<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Perc1<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      const auto IR  = static_cast<size_t>(Perc1InternalVariables::IR);
      const auto II  = static_cast<size_t>(Perc1InternalVariables::II);
      monitor_->set(Variable::ir, [this]
                    { return y_.getData()[IR]; });
      monitor_->set(Variable::ii, [this]
                    { return y_.getData()[II]; });
      monitor_->set(Variable::im, [this]
                    { return std::sqrt(y_.getData()[IR] * y_.getData()[IR] + y_.getData()[II] * y_.getData()[II]); });
      monitor_->set(Variable::p, [this]
                    { return bus_->Vr() * y_.getData()[IR] + bus_->Vi() * y_.getData()[II]; });
      monitor_->set(Variable::q, [this]
                    { return bus_->Vi() * y_.getData()[IR] - bus_->Vr() * y_.getData()[II]; });
      monitor_->set(Variable::ip, [this]
                    { return y_.getData()[static_cast<size_t>(Perc1InternalVariables::IP)]; });
      monitor_->set(Variable::iq, [this]
                    { return y_.getData()[static_cast<size_t>(Perc1InternalVariables::IQ)]; });
      monitor_->set(Variable::fracon, [this]
                    { return y_.getData()[static_cast<size_t>(Perc1InternalVariables::FRACON)]; });
    }
  } // namespace PhasorDynamics
} // namespace GridKit
