/**
 * @file ReecbImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT renewable energy electrical control.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Reecb/Reecb.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a renewable energy electrical control and its selections.
     *
     * Layout: 21 internal variables and inputs `vt`, `p`, `q`, `qext`,
     * `pfaref` and `pref`, each of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Reecb<scalar_type, index_type>::Reecb(const ModelDataT& data)
      : ComponentT(std::type_identity<Reecb>{}),
        Trv_(std::max(data.Trv, TIME_CONSTANT_MINIMUM)),
        Tp_(std::max(data.Tp, TIME_CONSTANT_MINIMUM)),
        Vref0_(data.Vref0),
        Vdip_(data.Vdip),
        Vup_(data.Vup),
        dbd1_(data.dbd1),
        dbd2_(data.dbd2),
        kqv_(data.kqv),
        Iql1_(data.Iql1),
        Iqh1_(data.Iqh1),
        Qmax_(data.Qmax),
        Qmin_(data.Qmin),
        Kqp_(data.Kqp),
        Kqi_(data.Kqi),
        Vmax_(data.Vmax),
        Vmin_(data.Vmin),
        Kvp_(data.Kvp),
        Kvi_(data.Kvi),
        Tiq_(std::max(data.Tiq, TIME_CONSTANT_MINIMUM)),
        Tpord_(std::max(data.Tpord, TIME_CONSTANT_MINIMUM)),
        dPmax_(data.dPmax),
        dPmin_(data.dPmin),
        Pmax_(data.Pmax),
        Pmin_(data.Pmin),
        Imax_(data.Imax),
        V_(data.V),
        mva_(data.mva)
    {
      if (data.PfFlag)
      {
        pf_on_ = ONE<RealT>;
      }
      if (data.QFlag)
      {
        q_on_ = ONE<RealT>;
      }
      if (data.QFlag && data.VFlag)
      {
        q_pi_on_ = ONE<RealT>;
      }
      if (data.QFlag && !data.VFlag)
      {
        v_ref_on_ = ONE<RealT>;
      }
      if (data.Pqflag)
      {
        pq_on_ = ONE<RealT>;
      }
      pf_off_   = ONE<RealT> - pf_on_;
      q_off_    = ONE<RealT> - q_on_;
      q_ref_on_ = ONE<RealT> - v_ref_on_;
      pq_off_   = ONE<RealT> - pq_on_;
    }

    /**
     * @brief Check the ratings, the voltage reference, the thresholds, the
     * deadband, the gains and the limits.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Reecb<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Reecb");
      checks.check(V_ > ZERO<RealT> && mva_ > ZERO<RealT>, "V and mva must be positive");
      checks.check(Vref0_ > ZERO<RealT>, "Vref0 must be positive: a value from the operating point is not supported");
      checks.check(Vdip_ < Vup_, "Vdip must be less than Vup");
      checks.check(dbd1_ <= ZERO<RealT> && ZERO<RealT> <= dbd2_, "dbd1 <= 0 <= dbd2 is required");
      checks.check(kqv_ >= ZERO<RealT>, "kqv must be nonnegative");
      checks.check(Iql1_ <= Iqh1_, "Iql1 must be less than or equal to Iqh1");
      checks.check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
      checks.check(Kqp_ >= ZERO<RealT> && Kqi_ >= ZERO<RealT>, "Kqp and Kqi must be nonnegative");
      checks.check(Vmin_ <= Vmax_, "Vmin must be less than or equal to Vmax");
      checks.check(Kvp_ >= ZERO<RealT> && Kvi_ >= ZERO<RealT>, "Kvp and Kvi must be nonnegative");
      checks.check(dPmin_ < ZERO<RealT> && ZERO<RealT> < dPmax_, "dPmin < 0 < dPmax is required");
      checks.check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
      checks.check(std::isfinite(Imax_) && Imax_ > ZERO<RealT>, "Imax must be finite and positive");
      return checks;
    }

    /**
     * @brief Rate limited between a negative and a positive bound.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] rate Rate without the limits.
     * @param[in] lower Negative bound.
     * @param[in] upper Positive bound.
     * @return The limited rate.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Reecb<scalar_type, index_type>::aslew(T rate, RealT lower, RealT upper)
    {
      return rate / (ONE<RealT> + Math::ramp(rate / upper - ONE<RealT>) + Math::ramp(rate / lower - ONE<RealT>));
    }

    /**
     * @brief Rate of a state that stays within a band around zero.
     *
     * The rate passes within the band, and only when it points inward at
     * either bound.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] state State.
     * @param[in] rate Rate without the limits.
     * @param[in] band Half width of the band.
     * @return The limited rate.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Reecb<scalar_type, index_type>::awband(T state, T rate, T band)
    {
      const T above_min = Math::above(state, -band);
      const T below_max = Math::below(state, band);
      return (above_min * below_max + (ONE<RealT> - below_max) * Math::sigmoid(-rate)
              + (ONE<RealT> - above_min) * Math::sigmoid(rate))
             * rate;
    }

    /**
     * @brief Square root where positive and zero elsewhere, with a finite slope at zero.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] x Argument.
     * @return The smooth square root.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Reecb<scalar_type, index_type>::sqrtramp(T x)
    {
      const RealT root_width = ONE<RealT> / Math::MU<RealT>;
      const RealT knee       = CURRENT_CIRCLE_KNEE;

      const T absolute    = std::abs(x);
      const T normalizer  = absolute + knee;
      const T scaled_x    = absolute / normalizer;
      const T scaled_knee = knee / normalizer;
      const T magnitude   = normalizer * std::sqrt(scaled_x * scaled_x + scaled_knee * scaled_knee);
      const T half_sum    = HALF<RealT> * magnitude + HALF<RealT> * absolute;
      const T conjugate   = QUARTER<RealT> * knee * (knee / half_sum);
      const T hinged      = HALF<RealT> * x + HALF<RealT> * absolute + conjugate;
      // Dividing by the root keeps a large finite circle from overflowing
      return hinged / std::sqrt(hinged + root_width * root_width);
    }

    /**
     * @brief Difference of the squares of a limit and of a current, in factored form.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] limit Current limit.
     * @param[in] high Current of the axis with priority.
     * @return `limit`^2 - `high`^2.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Reecb<scalar_type, index_type>::circleSquare(RealT limit, T high)
    {
      return ((HALF<RealT> * limit - HALF<RealT> * high) * (HALF<RealT> * limit + HALF<RealT> * high)) * FOUR<RealT>;
    }

    /**
     * @brief Write the residuals of the measurements, of the reactive and
     * active paths and of the current commands.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Reecb<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const RealT system_to_device = this->va_system_base_ / (MEGA<RealT> * mva_);

      const T vt   = w[VT_IN] / V_;
      const T pe   = w[P_IN] / (MEGA<RealT> * mva_);
      const T qgen = w[Q_IN] / (MEGA<RealT> * mva_);

      const T qref  = system_to_device * w[QEXT];
      const T pref  = system_to_device * w[PREF];
      const T sdip  = Math::inside(vt, Vdip_, Vup_);
      const T verr  = Math::deadband2(Vref0_ - y[VMEAS], dbd1_, dbd2_);
      const T q_pi  = Kqp_ * y[EQ] + y[XPIQ];
      const T v_pi  = Kvp_ * y[EPIV] + y[XPIV];
      const T fpord = (pref - y[PORD]) / Tpord_;
      // Selected before the factored square, so the inactive path never forms 0 * inf
      const T high  = pq_on_ * y[IPCMD] + pq_off_ * y[IQCMD];

      f[VMEAS]  = -yp[VMEAS] + (vt - y[VMEAS]) / Trv_;
      f[PMEAS]  = -yp[PMEAS] + (pe - y[PMEAS]) / Tp_;
      f[XPIQ]   = -yp[XPIQ] + q_pi_on_ * sdip * Math::antiwindup(q_pi, Kqi_ * y[EQ], Vmin_, Vmax_);
      f[XPIV]   = -yp[XPIV] + q_on_ * sdip * awband(v_pi, Kvi_ * y[EPIV], y[IQMAX]);
      f[QV]     = -yp[QV] + q_off_ * sdip * (y[QREF] / y[VSAFE] - y[QV]) / Tiq_;
      f[PORD]   = -yp[PORD] + sdip * Math::antiwindup(y[PORD], y[RPORD], Pmin_, Pmax_);
      f[VSAFE]  = -y[VSAFE] + Math::max(y[VMEAS], VMEAS_MINIMUM);
      f[SDIP]   = -y[SDIP] + sdip;
      f[IQV]    = -y[IQV] + Math::clamp(kqv_ * verr, Iql1_, Iqh1_);
      f[QREF]   = -y[QREF] + q_ref_on_ * (pf_on_ * y[PMEAS] * std::tan(w[PFAREF]) + pf_off_ * qref);
      f[EQ]     = -y[EQ] + Math::clamp(y[QREF], Qmin_, Qmax_) - qgen;
      f[VPIQ]   = -y[VPIQ] + Math::clamp(q_pi, Vmin_, Vmax_);
      f[EPIV]   = -y[EPIV] + q_pi_on_ * y[VPIQ] + v_ref_on_ * Math::clamp(w[QEXT], Vmin_, Vmax_) - q_on_ * y[VMEAS];
      f[RPORD]  = -y[RPORD] + aslew(fpord, dPmin_, dPmax_);
      f[ILCAP]  = -y[ILCAP] + sqrtramp(circleSquare(Imax_, high));
      f[IQMAX]  = -y[IQMAX] + pq_on_ * y[ILCAP] + pq_off_ * Imax_;
      f[IPMAX]  = -y[IPMAX] + pq_on_ * Imax_ + pq_off_ * y[ILCAP];
      f[IQBASE] = -y[IQBASE] + Math::clamp(v_pi, -y[IQMAX], y[IQMAX]);
      f[IQRAW]  = -y[IQRAW] + q_on_ * y[IQBASE] + q_off_ * y[QV] + y[IQV];
      f[IQCMD]  = -y[IQCMD] + Math::clamp(y[IQRAW], -y[IQMAX], y[IQMAX]);
      f[IPCMD]  = -y[IPCMD] + Math::clamp(y[PORD] / y[VSAFE], ZERO<RealT>, y[IPMAX]);
    }
  } // namespace EMT
} // namespace GridKit
