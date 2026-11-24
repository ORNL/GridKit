/**
 * @file ReecaImpl.hpp
 * @brief Equations of the WECC type A renewable electrical controller.
 */

#pragma once

#include <cmath>
#include <limits>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/CurrentCircle.hpp>
#include <GridKit/Model/EMT/Component/Controller/Reeca/Reeca.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    Reeca<ScalarT, IdxT>::Reeca(const ModelDataT& data)
      : ComponentT(std::type_identity<Reeca>{}), data_(data)
    {
      for (size_t k = 0; k < slope_q_.size(); ++k)
      {
        const auto slope = [k](const auto& v, const auto& i)
        {
          return v[k + 1] > v[k] ? (i[k + 1] - i[k]) / (v[k + 1] - v[k]) : RealT{0};
        };
        slope_q_[k] = slope(data.Vq, data.Iq);
        slope_p_[k] = slope(data.Vp, data.Ip);
      }
    }

    template <typename ScalarT, typename IdxT>
    Model::ConfigurationChecks Reeca<ScalarT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Reeca");
      const auto&                d = data_;
      checks.check(d.V > 0 && d.mva > 0, "V and mva must be positive");
      checks.check(d.Trv >= 0 && d.Tp >= 0 && d.Tiq >= 0 && d.Tpord > 0, "invalid control time constants");
      checks.check(d.Vref0 >= 0 && d.Vdip >= 0 && d.Vdip < d.Vup, "invalid reference or dip voltages");
      checks.check(d.Thld2 >= 0 && d.Imax > 0, "Thld2 must be nonnegative and Imax positive");
      checks.check(d.dbd1 <= 0 && d.dbd2 >= 0 && d.Iql1 <= d.Iqh1, "invalid injection limits or deadband");
      checks.check(d.dPmin < 0 && d.dPmax > 0 && d.Pmin <= d.Pmax, "invalid power or ramp limits");
      checks.check(d.Vmin <= d.Vmax && d.Qmin <= d.Qmax, "voltage and reactive limits must be ordered");
      const auto curve = [](const auto& v, const auto& i)
      {
        if (v == std::array<RealT, 4>{} && i == std::array<RealT, 4>{})
        {
          return true;
        }
        return v[0] >= 0 && v[0] < v[1] && v[1] < v[2] && v[2] < v[3] && std::ranges::min(i) >= 0;
      };
      checks.check(curve(d.Vq, d.Iq) && curve(d.Vp, d.Ip), "VDL tables must be empty or have increasing voltages and nonnegative limits");
      return checks;
    }

    template <typename ScalarT, typename IdxT>
    int Reeca<ScalarT, IdxT>::initialize()
    {
      const RealT voltage  = this->y_.getData()[VTERM];
      const RealT terminal = this->gather().first[VT_IN] / data_.V;
      const RealT sensed   = this->time() == 0 && terminal == 0 ? voltage : terminal;
      const bool  dip      = sensed < data_.Vdip || sensed > data_.Vup;
      if (this->time() == 0)
      {
        reference_voltage_ = data_.Vref0 > 0 ? data_.Vref0 : sensed;
        dip_               = dip;
        hold_p_ = hold_q_ = 0;
      }
      else
      {
        voltageDip(dip, this->time());
      }
      return 0;
    }

    template <typename ScalarT, typename IdxT>
    void Reeca<ScalarT, IdxT>::voltageDip(bool dip, RealT t)
    {
      if (dip == static_cast<bool>(dip_))
      {
        return;
      }
      dip_    = dip;
      hold_p_ = !dip && data_.Thld2 > 0;
      hold_q_ = !dip && data_.Thld != 0;
      if (!dip)
      {
        held_ipmax_ = this->y_.getData()[IPMAX];
        release_p_  = t + data_.Thld2;
        release_q_  = t + std::abs(data_.Thld);
      }
    }

    template <typename ScalarT, typename IdxT>
    int Reeca<ScalarT, IdxT>::evaluateRootsAt(RealT t, const ScalarT* y, RealT* roots)
    {
      roots[0] = y[VTERM] - data_.Vdip;
      roots[1] = y[VTERM] - data_.Vup;
      roots[2] = hold_p_ != 0 ? t - release_p_ : RealT{1};
      roots[3] = hold_q_ != 0 ? t - release_q_ : RealT{1};
      return 0;
    }

    template <typename ScalarT, typename IdxT>
    bool Reeca<ScalarT, IdxT>::applyRootEvents(RealT t, std::span<const int> crossed)
    {
      const bool was_dipped = dip_ != 0;
      if (crossed[2])
      {
        hold_p_ = 0;
      }
      if (crossed[3])
      {
        hold_q_ = 0;
      }
      if (crossed[0] || crossed[1])
      {
        voltageDip(crossed[0] < 0 || crossed[1] > 0, t);
      }
      else
      {
        // Co-located controllers can reach the same threshold within IDA's root-time roundoff.
        const RealT v     = this->y_.getData()[VTERM];
        const RealT dv    = this->yp_.getData()[VTERM];
        const RealT round = 100 * std::numeric_limits<RealT>::epsilon() * (1 + std::abs(v) + std::abs(t * dv));
        if (dv != 0 && std::abs(v - data_.Vdip) <= round)
        {
          voltageDip(dv < 0, t);
        }
        else if (dv != 0 && std::abs(v - data_.Vup) <= round)
        {
          voltageDip(dv > 0, t);
        }
      }
      return was_dipped != (dip_ != 0) || crossed[2] || crossed[3];
    }

    template <typename ScalarT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] void
    Reeca<ScalarT, IdxT>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto& d      = data_;
      const RealT base   = MEGA<RealT> * d.mva;
      const RealT ratio  = this->va_system_base_ / base;
      const T     speed  = ONE<RealT> + d.Pflag * (w[OMEGA_IN] - ONE<RealT>);
      const T     pref   = ratio * w[PREF_IN] * speed;
      const T     qref   = d.PfFlag * y[PMEAS] * std::tan(w[PFAREF_IN]) + !d.PfFlag * ratio * w[QEXT_IN];
      const T     wanted = Math::clamp(qref, d.Qmin, d.Qmax);
      const T     vpi    = d.Kqp * y[EQ] + y[XPIQ];
      const T     iqpi   = d.Kvp * y[EPIV] + y[XPIV];
      const auto  limit  = [&](const auto& v, const auto& i, const auto& slope)
      {
        T current = i.front() + (v.back() == 0) * d.Imax;
#pragma clang loop unroll(full)
        for (size_t k = 0; k < 3; ++k)
        {
          current += slope[k] * (Math::ramp(y[VMEAS] - v[k]) - Math::ramp(y[VMEAS] - v[k + 1]));
        }
        return current - std::fdim(current, d.Imax);
      };
      const T iq_table = limit(d.Vq, d.Iq, slope_q_);
      const T ip_table = limit(d.Vp, d.Ip, slope_p_);
      const T ipmax    = ip_table - !d.Pqflag * std::fdim(ip_table, y[ILCAP]);
      const T iq       = d.QFlag * Math::clamp(iqpi, -y[IQMAX], y[IQMAX]) + !d.QFlag * y[QV]
                   + dip_ * y[IQV] + hold_q_ * ((d.Thld < 0) * y[IQV] + (d.Thld > 0) * d.Iqfrz);
      const RealT     run = ONE<RealT> - dip_;
      constexpr RealT tracking_time{1.0e-3};
      const T         tracking_rate = (pref - y[PREF_RATE]) / tracking_time;

      f[VSAFE] = Math::max(y[VMEAS], RealT{0.01}) - y[VSAFE];
      f[EQ]    = wanted - w[Q_IN] / base - y[EQ];
      f[EPIV]  = d.VFlag * Math::clamp(vpi, d.Vmin, d.Vmax)
                + !d.VFlag * (w[QEXT_IN] + d.Vref1) - y[VMEAS] - y[EPIV];
      f[ILCAP]     = currentCircle(d.Imax, d.Pqflag * y[IPCMD] + !d.Pqflag * y[IQCMD]) - y[ILCAP];
      f[IQMAX]     = iq_table - d.Pqflag * std::fdim(iq_table, y[ILCAP]) - y[IQMAX];
      f[VTERM]     = w[VT_IN] / d.V - y[VTERM];
      f[VMEAS]     = y[VTERM] - y[VMEAS] - d.Trv * yp[VMEAS];
      f[PMEAS]     = w[P_IN] / base - y[PMEAS] - d.Tp * yp[PMEAS];
      f[XPIQ]      = (d.QFlag && d.VFlag) * run * Math::antiwindup(vpi, d.Kqi * y[EQ], d.Vmin, d.Vmax) - yp[XPIQ];
      f[XPIV]      = d.QFlag * run * Math::antiwindup(iqpi, d.Kvi * y[EPIV], -y[IQMAX], y[IQMAX]) - yp[XPIV];
      f[QV]        = (ONE<RealT> - (d.Tiq > 0) * dip_) * (wanted / y[VSAFE] - y[QV]) - d.Tiq * yp[QV];
      f[PREF_RATE] = std::clamp(tracking_rate, T{d.dPmin}, T{d.dPmax}) - yp[PREF_RATE];
      f[PORD]      = run * Math::antiwindup(y[PORD], (y[PREF_RATE] - y[PORD]) / d.Tpord, d.Pmin, d.Pmax) - yp[PORD];
      f[IQV]       = Math::clamp(d.kqv * Math::deadband2(reference_voltage_ - y[VMEAS], d.dbd1, d.dbd2), d.Iql1, d.Iqh1) - y[IQV];
      f[IPMAX]     = (ONE<RealT> - hold_p_) * ipmax + hold_p_ * held_ipmax_ - y[IPMAX];
      f[IQCMD]     = Math::clamp(iq, -y[IQMAX], y[IQMAX]) - y[IQCMD];
      f[IPCMD]     = Math::clamp(y[PORD] / y[VSAFE], ZERO<RealT>, y[IPMAX]) - y[IPCMD];
    }
  } // namespace EMT
} // namespace GridKit
