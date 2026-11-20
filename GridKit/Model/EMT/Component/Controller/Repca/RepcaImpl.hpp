/**
 * @file RepcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT renewable energy plant control.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Repca/Repca.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a renewable energy plant control and its selections.
     *
     * Layout: 21 internal variables and inputs `v`, `p`, `q`, `omega`, `vref`,
     * `pref`, `qref` and `freqref`, each of size 1. `verify` reports invalid
     * data.
     */
    template <typename scalar_type, typename index_type>
    Repca<scalar_type, index_type>::Repca(const ModelDataT& data)
      : ComponentT(std::type_identity<Repca>{}),
        mva_(data.mva),
        Tfltr_(std::max(data.Tfltr, TIME_CONSTANT_MINIMUM)),
        Vfrz_(data.Vfrz),
        Rc_(data.Rc),
        Xc_(data.Xc),
        Kc_(data.Kc),
        dbdlow_(data.dbdlow),
        dbdupper_(data.dbdupper),
        emax_(data.emax),
        emin_(data.emin),
        Kp_(data.Kp),
        Ki_(data.Ki),
        Qmax_(data.Qmax),
        Qmin_(data.Qmin),
        Tft_(data.Tft),
        Tfv_(std::max(data.Tfv, TIME_CONSTANT_MINIMUM)),
        Tp_(std::max(data.Tp, TIME_CONSTANT_MINIMUM)),
        fdbd1_(data.fdbd1),
        fdbd2_(data.fdbd2),
        Ddn_(data.Ddn),
        Dup_(data.Dup),
        femax_(data.femax),
        femin_(data.femin),
        Kpg_(data.Kpg),
        Kig_(data.Kig),
        Pmax_(data.Pmax),
        Pmin_(data.Pmin),
        Tlag_(std::max(data.Tlag, TIME_CONSTANT_MINIMUM)),
        V_(data.V)
    {
      if (data.VcompFlag)
      {
        vcomp_on_ = ONE<RealT>;
      }
      if (data.RefFlag)
      {
        ref_on_ = ONE<RealT>;
      }
      if (data.Freqflag)
      {
        freq_on_ = ONE<RealT>;
      }
      vcomp_off_ = ONE<RealT> - vcomp_on_;
      ref_off_   = ONE<RealT> - ref_on_;

      // With a proportional gain the integrator tracks the limited output
      q_integral_gain_ = Ki_;
      if (Kp_ > ZERO<RealT>)
      {
        q_feedback_gain_ = Ki_ / Kp_;
        q_integral_gain_ = ZERO<RealT>;
      }
      p_integral_gain_ = Kig_;
      if (Kpg_ > ZERO<RealT>)
      {
        p_feedback_gain_ = Kig_ / Kpg_;
        p_integral_gain_ = ZERO<RealT>;
      }
    }

    /**
     * @brief Check the ratings, the gains, the deadbands and the limits.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Repca<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Repca");
      checks.check(V_ > ZERO<RealT> && mva_ > ZERO<RealT>, "V and mva must be positive");
      checks.check(Kp_ >= ZERO<RealT>, "Kp must be nonnegative");
      checks.check(Kpg_ >= ZERO<RealT>, "Kpg must be nonnegative");
      checks.check(dbdlow_ <= ZERO<RealT> && ZERO<RealT> <= dbdupper_, "dbdlow <= 0 <= dbdupper is required");
      checks.check(emin_ <= ZERO<RealT> && ZERO<RealT> <= emax_, "emin <= 0 <= emax is required");
      checks.check(Qmin_ <= Qmax_, "Qmin must be less than or equal to Qmax");
      checks.check(fdbd1_ <= ZERO<RealT> && ZERO<RealT> <= fdbd2_, "fdbd1 <= 0 <= fdbd2 is required");
      checks.check(Ddn_ >= ZERO<RealT>, "Ddn must be nonnegative");
      checks.check(Dup_ >= ZERO<RealT>, "Dup must be nonnegative");
      checks.check(femin_ <= ZERO<RealT> && ZERO<RealT> <= femax_, "femin <= 0 <= femax is required");
      checks.check(Pmin_ <= Pmax_, "Pmin must be less than or equal to Pmax");
      return checks;
    }

    /**
     * @brief Frequency response with one gain for each sign of the error.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] error Frequency error.
     * @param[in] down Gain for a negative error.
     * @param[in] up Gain for a positive error.
     * @return The response.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Repca<scalar_type, index_type>::droop(T error, RealT down, RealT up)
    {
      return error * (down + (up - down) * Math::sigmoid(error));
    }

    /**
     * @brief Write the residuals of the measurements, of the reactive and
     * active loops and of the power commands.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Repca<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const RealT system_to_device = this->va_system_base_ / (MEGA<RealT> * mva_);
      const RealT omega_base       = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;

      const T vmag = w[V_IN] / V_;
      const T p    = w[P_IN] / (MEGA<RealT> * mva_);
      const T q    = w[Q_IN] / (MEGA<RealT> * mva_);

      // Squared magnitude of the voltage less the drop across Rc + j Xc, whose current has the square (p^2 + q^2) / v^2
      const T i2    = (p * p + q * q) / Math::max(vmag * vmag, V_MINIMUM * V_MINIMUM);
      const T vldc2 = vmag * vmag - TWO<RealT> * (Rc_ * p + Xc_ * q) + (Rc_ * Rc_ + Xc_ * Xc_) * i2;
      T       vldc  = ZERO<RealT> * vldc2;
      if (vldc2 > ZERO<RealT>)
      {
        vldc = std::sqrt(vldc2);
      }

      const T sfrz      = Math::above(vmag, Vfrz_);
      const T qext      = system_to_device * y[QEXT];
      const T pext      = system_to_device * y[PEXT];
      const T pfreq     = droop(y[EF], Ddn_, Dup_);
      const T q_pi_rate = q_feedback_gain_ * (y[QPI] - Math::clamp(y[XQPI], Qmin_, Qmax_))
                          + Math::antiwindup(y[XQPI], q_integral_gain_ * y[ERQLIM], Qmin_, Qmax_);
      const T p_pi_rate = p_feedback_gain_ * (y[PPI] - Math::clamp(y[XPPI], Pmin_, Pmax_))
                          + Math::antiwindup(y[XPPI], p_integral_gain_ * y[EPLIM], Pmin_, Pmax_);

      f[VMEAS]  = -yp[VMEAS] + (y[VCTRL] - y[VMEAS]) / Tfltr_;
      f[QMEAS]  = -yp[QMEAS] + (q - y[QMEAS]) / Tfltr_;
      f[XQPI]   = -yp[XQPI] + sfrz * q_pi_rate;
      f[XQLAG]  = -yp[XQLAG] + (y[QPI] - y[XQLAG]) / Tfv_;
      f[PMEAS]  = -yp[PMEAS] + (p - y[PMEAS]) / Tp_;
      f[XPPI]   = -yp[XPPI] + p_pi_rate;
      f[PREF]   = -yp[PREF] + (y[PPI] - y[PREF]) / Tlag_;
      f[VLDC]   = -y[VLDC] + vldc;
      f[VDROOP] = -y[VDROOP] + vmag + Kc_ * q;
      f[VCTRL]  = -y[VCTRL] + vcomp_on_ * y[VLDC] + vcomp_off_ * y[VDROOP];
      f[SFRZ]   = -y[SFRZ] + sfrz;
      f[ERQ]    = -y[ERQ] + ref_on_ * (w[VREF] - y[VMEAS]) + ref_off_ * (system_to_device * w[QREF] - y[QMEAS]);
      f[ERQDB]  = -y[ERQDB] + Math::deadband2(y[ERQ], dbdlow_, dbdupper_);
      f[ERQLIM] = -y[ERQLIM] + Math::clamp(y[ERQDB], emin_, emax_);
      f[QPI]    = -y[QPI] + Math::clamp(Kp_ * y[ERQLIM] + y[XQPI], Qmin_, Qmax_);
      f[QEXT]   = -Tfv_ * (qext - y[XQLAG]) + Tft_ * (y[QPI] - y[XQLAG]);
      f[EF]     = -y[EF] + Math::deadband2(w[FREQREF] - w[OMEGA] / omega_base, fdbd1_, fdbd2_);
      f[EP]     = -y[EP] + system_to_device * w[PREF_IN] - y[PMEAS] + pfreq;
      f[EPLIM]  = -y[EPLIM] + Math::clamp(y[EP], femin_, femax_);
      f[PPI]    = -y[PPI] + Math::clamp(Kpg_ * y[EPLIM] + y[XPPI], Pmin_, Pmax_);
      f[PEXT]   = -pext + freq_on_ * y[PREF];
    }
  } // namespace EMT
} // namespace GridKit
