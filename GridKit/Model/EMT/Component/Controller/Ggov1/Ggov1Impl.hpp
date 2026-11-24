/**
 * @file Ggov1Impl.hpp
 * @brief Equations of the general turbine governor.
 */

#pragma once

#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ggov1/Ggov1.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    Ggov1<ScalarT, IdxT>::Ggov1(const ModelDataT& data)
      : ComponentT(std::type_identity<Ggov1>{}), data_(data)
    {
      const RealT rating = data.Trate > 0 ? data.Trate * MEGA<RealT> : data.S;
      if (rating > 0 && data.S > 0)
      {
        inv_rating_ = ONE<RealT> / rating;
        to_machine_ = rating / data.S;
      }
      governor_tracking_ = data.Kpgov > 0 ? data.Kigov / data.Kpgov : 0;
      load_tracking_     = data.Kpload > 0 ? data.Kiload / data.Kpload : 0;
      temperature_gain_  = data.Tsb > 0 ? data.Tsa / data.Tsb : 0;
      turbine_gain_      = data.Tb > 0 ? data.Tc / data.Tb : 0;
      acceleration_gain_ = data.Ta > 0 ? data.Ka / data.Ta : 0;
      derivative_gain_   = data.Tdgov > 0 ? data.Kdgov / data.Tdgov : 0;
    }

    template <typename ScalarT, typename IdxT>
    Model::ConfigurationChecks Ggov1<ScalarT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Ggov1");
      const auto&                d = data_;
      checks.check(d.S > 0 && d.Trate >= 0, "S must be positive and Trate nonnegative");
      checks.check(-2 <= d.Rselect && d.Rselect <= 1 && d.R >= 0, "invalid droop selection or gain");
      checks.check(d.Kpgov >= 0 && d.Kigov >= 0 && d.Kpload >= 0 && d.Kiload >= 0, "PI gains must be nonnegative");
      checks.check(d.Tpelec >= 0 && d.Tdgov >= 0 && (d.Tdgov > 0 || d.Kdgov == 0), "invalid measurement or derivative time");
      checks.check(d.Tact > 0 && d.Kturb > 0 && d.Tfload >= 0, "invalid actuator, turbine or load limiter");
      checks.check(d.Tb >= 0 && d.Tc >= 0 && (d.Tb > 0 || d.Tc == 0), "Tc requires a positive Tb");
      checks.check(d.Tsb >= 0 && d.Tsa >= 0 && (d.Tsb > 0 || d.Tsa == 0), "Tsa requires a positive Tsb");
      checks.check(d.Teng == 0 && d.Tfsr > 0, "GGOV1 requires zero Teng and positive fuel memory Tfsr");
      checks.check(d.Ta >= 0 && (d.Ka == 0 || d.Ta > 0), "acceleration control requires a positive Ta");
      checks.check(d.Minerr <= d.Maxerr && d.Vmin <= d.Vmax && d.Rclose < 0 && d.Ropen > 0, "invalid governor limits");
      return checks;
    }

    template <typename ScalarT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] void
    Ggov1<ScalarT, IdxT>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto& d        = data_;
      const T     speed    = w[OMEGA_IN] - ONE<RealT>;
      const T     feedback = (d.Rselect == 1) * y[PE] + (d.Rselect == -1) * y[VALVE] + (d.Rselect == -2) * y[FSR];
      f[ERROR]             = Math::clamp(Math::deadband2(w[PREF_IN] + y[XMW] - speed - d.R * feedback, -d.Db, d.Db), d.Minerr, d.Maxerr) - y[ERROR];
      const T derivative   = derivative_gain_ * (y[ERROR] - y[XD]);
      f[FSRN]              = d.Kpgov * y[ERROR] + y[XI] + derivative - y[FSRN];
      const T fuel         = y[VALVE] * (ONE<RealT> + d.Flag * speed);
      const T exhaust      = fuel * std::exp(std::min(d.Dm, RealT{0}) * std::log(Math::max(w[OMEGA_IN], RealT{0.01})));
      const T sensed       = y[TEMP] + temperature_gain_ * (exhaust - y[TEMP]);
      const T load_error   = (d.Ldref + w[PAUX_IN]) / d.Kturb + d.Wfnl - y[LOAD];
      f[FSRT]              = d.Kpload * load_error + y[XL] - y[FSRT];
      const T turbine      = y[TURB] + turbine_gain_ * (fuel - y[TURB]);

      f[PE]               = w[P_IN] * inv_rating_ - y[PE] - d.Tpelec * yp[PE];
      f[XD]               = y[ERROR] - y[XD] - d.Tdgov * yp[XD];
      f[XI]               = governor_tracking_ * (y[FSR] - y[XI]) + (d.Kpgov == 0) * d.Kigov * y[ERROR] - yp[XI];
      f[VALVE]            = Math::clamp((y[FSR] - y[VALVE]) / d.Tact, d.Rclose, d.Ropen) - yp[VALVE];
      f[TURB]             = fuel - y[TURB] - d.Tb * yp[TURB];
      f[TEMP]             = exhaust - y[TEMP] - d.Tsb * yp[TEMP];
      f[LOAD]             = sensed - y[LOAD] - d.Tfload * yp[LOAD];
      f[XL]               = load_tracking_ * (y[FSR] - y[XL]) + (d.Kpload == 0) * d.Kiload * load_error - yp[XL];
      f[XMW]              = Math::antiwindup(y[XMW], d.Kimw * (w[PMWSET_IN] - y[PE]), -1.1 * d.R, 1.1 * d.R) - yp[XMW];
      f[SPEED]            = speed - y[SPEED] - d.Ta * yp[SPEED];
      // Exact fuel selection preserves the operating point at coincident orders and valve limits.
      const auto positive = [](T x)
      { return HALF<RealT> * (x + std::abs(x)); };
      const auto minimum = [&](T a, T b)
      { return a - positive(a - b); };
      const auto limited = [&](T x)
      { return d.Vmin + positive(x - d.Vmin) - positive(x - d.Vmax); };
      f[DEMAND]     = limited(minimum(y[FSRN], minimum(y[FSRT], d.Fsrtmax))) - y[DEMAND];
      // The acceleration controller compares demand with the previous fuel command.
      const T accel = d.Ka * d.Aset - acceleration_gain_ * (speed - y[SPEED]);
      const T order = y[DEMAND] - (d.Ka > 0) * positive(y[DEMAND] - y[MEMORY] - d.Tfsr * accel);
      f[FSR]        = limited(order) - y[FSR];
      f[MEMORY]     = y[FSR] - y[MEMORY] - d.Tfsr * yp[MEMORY];
      f[PM]         = to_machine_ * (d.Kturb * (turbine - d.Wfnl) - std::max(d.Dm, RealT{0}) * speed) - y[PM];
    }
  } // namespace EMT
} // namespace GridKit
