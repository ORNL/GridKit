/**
 * @file Esdc1aImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT IEEE type DC1A excitation system.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esdc1a/Esdc1a.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct an IEEE type DC1A excitation system, its selections and its saturation.
     *
     * Layout: 11 internal variables of size 1, and inputs `vt`, `omega`,
     * `vref`, `vs` and `vuel` of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Esdc1a<scalar_type, index_type>::Esdc1a(const ModelDataT& data)
      : ComponentT(std::type_identity<Esdc1a>{}),
        V_(data.V),
        Tr_(std::max(data.Tr, TIME_CONSTANT_MINIMUM)),
        Ka_(data.Ka),
        Ta_(std::max(data.Ta, TIME_CONSTANT_MINIMUM)),
        Tb_(std::max(data.Tb, TIME_CONSTANT_MINIMUM)),
        Tc_(data.Tc),
        Vrmax_(data.Vrmax),
        Vrmin_(data.Vrmin),
        Ke_(data.Ke),
        Te_(std::max(data.Te, TIME_CONSTANT_MINIMUM)),
        Kf_(data.Kf),
        Tf1_(std::max(data.Tf1, TIME_CONSTANT_MINIMUM)),
        E1_(data.E1),
        Se1_(data.Se1),
        E2_(data.E2),
        Se2_(data.Se2),
        UEL_(data.UEL)
    {
      if (data.Spdmlt)
      {
        spd_on_ = ONE<RealT>;
      }
      if (UEL_ >= SUMMING_UEL)
      {
        uel_on_ = ONE<RealT>;
      }
      if (data.exclim)
      {
        lim_on_ = ONE<RealT>;
      }

      // Invalid or absent saturation data leaves SA and SB zero
      const bool ordered = (E2_ > E1_ && Se2_ > Se1_) || (E2_ < E1_ && Se2_ < Se1_);
      if (!ordered || E1_ <= ZERO<RealT> || E2_ <= ZERO<RealT> || Se1_ < ZERO<RealT> || Se2_ < ZERO<RealT>)
      {
        return;
      }
      if (Se1_ == ZERO<RealT>)
      {
        SA_ = E1_;
        SB_ = Se2_ * E2_ / ((E2_ - E1_) * (E2_ - E1_));
        return;
      }
      if (Se2_ == ZERO<RealT>)
      {
        SA_ = E2_;
        SB_ = Se1_ * E1_ / ((E1_ - E2_) * (E1_ - E2_));
        return;
      }
      const RealT C = std::sqrt(Se2_ * E2_ / (Se1_ * E1_));
      SA_           = (C * E1_ - E2_) / (C - ONE<RealT>);
      SB_           = Se1_ * E1_ / ((E1_ - SA_) * (E1_ - SA_));
    }

    /**
     * @brief Check the rated voltage, the gains, the regulator limits, the
     * limiter selection and the saturation data.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Esdc1a<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Esdc1a");
      checks.check(V_ > ZERO<RealT>, "V must be positive");
      checks.check(Ka_ > ZERO<RealT>, "Ka must be positive");
      checks.check(Ke_ != ZERO<RealT>, "Ke must be nonzero: a value from the operating point is not supported");
      checks.check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
      checks.check(UEL_ <= MAX_UEL, "UEL must be 0, 1, 2 or 3");
      if (Se1_ != ZERO<RealT> || Se2_ != ZERO<RealT>)
      {
        checks.check(E1_ > ZERO<RealT> && E2_ > ZERO<RealT>, "E1 and E2 must be positive with saturation");
        checks.check(Se1_ >= ZERO<RealT> && Se2_ >= ZERO<RealT>, "Se1 and Se2 must be nonnegative with saturation");
        checks.check((E2_ > E1_ && Se2_ > Se1_) || (E2_ < E1_ && Se2_ < Se1_),
                     "E1, E2 and Se1, Se2 must be ordered consistently");
      }
      return checks;
    }

    /**
     * @brief Rate of a state that stays above a lower limit.
     *
     * The rate passes above the limit, and only when positive at it.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] x State.
     * @param[in] f Rate without the limit.
     * @param[in] xmin Lower limit.
     * @return The limited rate.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Esdc1a<scalar_type, index_type>::awmin(T x, T f, RealT xmin)
    {
      const T above = Math::above(x, xmin);
      return (above + (ONE<RealT> - above) * Math::sigmoid(f)) * f;
    }

    /**
     * @brief Write the residuals of the voltage sensing, of the lead-lag, of
     * the regulator, of the exciter and of the stabilizing feedback.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Esdc1a<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T vt           = w[VT] / V_;
      const T omega        = w[OMEGA] - ONE<RealT>;
      const T rate         = (y[VR] - y[VFE]) / Te_;
      const T limited_rate = (ONE<RealT> - lim_on_) * rate + lim_on_ * awmin(y[EFDP], rate, ZERO<RealT>);

      f[EFDP] = -yp[EFDP] + limited_rate;
      f[VC]   = -yp[VC] + (vt - y[VC]) / Tr_;
      f[VR]   = -yp[VR] + Math::antiwindup(y[VR], -y[VR] + Ka_ * y[VHV], Vrmin_, Vrmax_) / Ta_;
      f[VF]   = -yp[VF] + (-y[VF] + Kf_ * limited_rate) / Tf1_;
      f[XLL]  = -yp[XLL] + (y[EV] - y[XLL]) / Tb_;
      f[EV]   = -y[EV] + w[VREF] + w[VS] + uel_on_ * w[VUEL] - y[VC] - y[VF];
      f[VLL]  = -y[VLL] + y[XLL] + (Tc_ / Tb_) * (y[EV] - y[XLL]);
      f[VHV]  = -y[VHV] + uel_on_ * y[VLL] + (ONE<RealT> - uel_on_) * Math::max(y[VLL], w[VUEL]);
      f[SE]   = -y[SE] + SB_ * Math::qramp(y[EFDP] - SA_);
      f[VFE]  = -y[VFE] + Ke_ * y[EFDP] + y[SE];
      f[EFD]  = -y[EFD] + (ONE<RealT> + spd_on_ * omega) * y[EFDP];
    }
  } // namespace EMT
} // namespace GridKit
