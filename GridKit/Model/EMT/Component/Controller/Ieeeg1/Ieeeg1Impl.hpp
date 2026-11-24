/**
 * @file Ieeeg1Impl.hpp
 * @brief Equations of the IEEE type 1 steam turbine governor.
 */

#pragma once

#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeeg1/Ieeeg1.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    Ieeeg1<ScalarT, IdxT>::Ieeeg1(const ModelDataT& data)
      : ComponentT(std::type_identity<Ieeeg1>{}), data_(data)
    {
      to_machine_     = data.S > 0 && data.Trate > 0 ? MEGA<RealT> * data.Trate / data.S : ONE<RealT>;
      lead_gain_      = data.T1 > 0 ? data.T2 / data.T1 : ZERO<RealT>;
      straight_valve_ = data.Gv == decltype(data.Gv){} && data.Pgv == decltype(data.Pgv){};
      for (size_t k = 0; k < slope_.size(); ++k)
      {
        if (data.Gv[k + 1] > data.Gv[k])
        {
          slope_[k] = (data.Pgv[k + 1] - data.Pgv[k]) / (data.Gv[k + 1] - data.Gv[k]);
        }
      }
    }

    template <typename ScalarT, typename IdxT>
    Model::ConfigurationChecks Ieeeg1<ScalarT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Ieeeg1");
      const auto&                d = data_;
      checks.check(d.S > 0 && d.Trate >= 0, "S must be positive and Trate nonnegative");
      checks.check(d.T1 >= 0 && d.T2 >= 0 && (d.T1 > 0 || d.T2 == 0), "T2 requires a positive T1");
      checks.check(d.T3 > 0 && d.T4 >= 0 && d.T5 >= 0 && d.T6 >= 0 && d.T7 >= 0, "invalid servo or stage time constant");
      checks.check(d.Uc < 0 && d.Uo > 0 && d.Pmin <= d.Pmax, "invalid valve limits");
      checks.check(d.db1 >= 0 && d.Eps == 0 && d.db2 == 0, "db1 must be nonnegative; hysteresis and backlash require zero Eps and db2");
      if (d.Gv != decltype(d.Gv){} || d.Pgv != decltype(d.Pgv){})
      {
        for (size_t k = 1; k < d.Gv.size(); ++k)
        {
          checks.check(d.Gv[k] > d.Gv[k - 1] && d.Pgv[k] >= d.Pgv[k - 1], "valve curve must be monotone with rising Gv");
        }
      }
      return checks;
    }

    template <typename ScalarT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] void
    Ieeeg1<ScalarT, IdxT>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto& d         = data_;
      const T     deviation = w[OMEGA_IN] - ONE<RealT>;
      const T     speed     = d.K * Math::deadband1(deviation, -d.db1, d.db1);
      const T     lead      = y[SPEED] + lead_gain_ * (speed - y[SPEED]);
      const T     rate      = (w[PREF_IN] + w[PAUX_IN] - lead - y[VALVE]) / d.T3;
      T           flow      = straight_valve_ * y[VALVE] + d.Pgv.front();
#pragma clang loop unroll(full)
      for (size_t k = 0; k < slope_.size(); ++k)
      {
        flow += slope_[k] * (Math::ramp(y[VALVE] - d.Gv[k]) - Math::ramp(y[VALVE] - d.Gv[k + 1]));
      }
      f[SPEED]     = speed - y[SPEED] - d.T1 * yp[SPEED];
      f[VALVE]     = Math::boundedAntiwindup(y[VALVE], Math::clamp(rate, d.Uc, d.Uo), d.Pmin, d.Pmax) - yp[VALVE];
      f[BOWL]      = flow - y[BOWL] - d.T4 * yp[BOWL];
      f[REHEAT]    = y[BOWL] - y[REHEAT] - d.T5 * yp[REHEAT];
      f[CROSSOVER] = y[REHEAT] - y[CROSSOVER] - d.T6 * yp[CROSSOVER];
      f[REHEAT2]   = y[CROSSOVER] - y[REHEAT2] - d.T7 * yp[REHEAT2];
      f[PM]        = to_machine_ * (d.K1 * y[BOWL] + d.K3 * y[REHEAT] + d.K5 * y[CROSSOVER] + d.K7 * y[REHEAT2]) - y[PM];
      f[PM_LP]     = to_machine_ * (d.K2 * y[BOWL] + d.K4 * y[REHEAT] + d.K6 * y[CROSSOVER] + d.K8 * y[REHEAT2]) - y[PM_LP];
    }
  } // namespace EMT
} // namespace GridKit
