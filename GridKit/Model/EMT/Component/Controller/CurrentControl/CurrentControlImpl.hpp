/**
 * @file CurrentControlImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT current control.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/CurrentControl/CurrentControl.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a current control and its peak current.
     *
     * Layout: internal variables `xdq` and `mdq` of size 2, inputs `vpos`,
     * `vneg`, `omega`, `ip` and `iq` of size 1, and inputs `vdq`, `idq` and
     * `axis` of size 2. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    CurrentControl<scalar_type, index_type>::CurrentControl(const ModelDataT& data)
      : ComponentT(std::type_identity<CurrentControl>{}),
        V_(data.V),
        mva_(data.mva),
        Kp_(data.Kp),
        Ki_(data.Ki),
        L_(data.L),
        C_(data.C),
        Mmax_(data.Mmax)
    {
      if (V_ > ZERO<RealT>)
      {
        current_peak_ = std::sqrt(TWO<RealT> / THREE<RealT>) * MEGA<RealT> * mva_ / V_;
      }
    }

    /**
     * @brief Check the ratings, the gains, the inductance, the capacitance
     *        and the modulation limit.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks CurrentControl<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("CurrentControl");
      checks.check(V_ > ZERO<RealT> && mva_ > ZERO<RealT> && std::isfinite(current_peak_), "V and mva must be positive");
      checks.check(Kp_ > ZERO<RealT>, "Kp must be positive");
      checks.check(Ki_ > ZERO<RealT>, "Ki must be positive");
      checks.check(L_ > ZERO<RealT>, "L must be positive");
      checks.check(std::isfinite(C_) && C_ >= ZERO<RealT>, "C must be finite and nonnegative");
      checks.check(Mmax_ > ZERO<RealT>, "Mmax must be positive");
      return checks;
    }

    /**
     * @brief Set the integrator tolerances to the tolerance of the positive
     *        DC terminal voltage, and the modulation tolerances to `rel_tol`.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre What input `vpos` reads has set its tolerance.
     */
    template <typename scalar_type, typename index_type>
    int CurrentControl<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const RealT v_tol = this->externalTolerance(this->inputOffset(CurrentControlInputs::vpos));

      std::ranges::fill(this->absoluteTolerance(CurrentControlInternalVariables::xdq), v_tol);
      std::ranges::fill(this->absoluteTolerance(CurrentControlInternalVariables::mdq), rel_tol);
      return 0;
    }

    /**
     * @brief Write the residuals of the integrator states and of the
     *        modulation signals.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    CurrentControl<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const RealT scale = TWO<RealT> / Mmax_;

      const T id_cmd = current_peak_ * (w[AXIS_D] * w[IP_CMD] + w[AXIS_Q] * w[IQ_CMD]);
      const T iq_cmd = current_peak_ * (w[AXIS_Q] * w[IP_CMD] - w[AXIS_D] * w[IQ_CMD]);
      const T vdc    = w[VPOS] - w[VNEG];
      const T ed     = id_cmd - w[OMEGA] * C_ * w[VQ] - w[ID];
      const T eq     = iq_cmd + w[OMEGA] * C_ * w[VD] - w[IQ];
      const T ud     = w[VD] - w[OMEGA] * L_ * w[IQ] + Kp_ * ed + y[XD];
      const T uq     = w[VQ] + w[OMEGA] * L_ * w[ID] + Kp_ * eq + y[XQ];

      const T vmod   = std::sqrt(Math::max(vdc * vdc, scale * scale * (ud * ud + uq * uq)));
      const T ud_lim = vdc / vmod * ud;
      const T uq_lim = vdc / vmod * uq;

      f[XD] = -yp[XD] + Ki_ * (ed + (ud_lim - ud) / Kp_);
      f[XQ] = -yp[XQ] + Ki_ * (eq + (uq_lim - uq) / Kp_);
      f[MD] = -y[MD] + TWO<RealT> / vmod * ud;
      f[MQ] = -y[MQ] + TWO<RealT> / vmod * uq;
    }
  } // namespace EMT
} // namespace GridKit
