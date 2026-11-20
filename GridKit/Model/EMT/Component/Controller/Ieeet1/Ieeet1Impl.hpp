/**
 * @file Ieeet1Impl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT IEEE type 1 excitation system.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeet1/Ieeet1.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct an IEEE type 1 excitation system and its saturation.
     *
     * Layout: 9 internal variables of size 1, and inputs `vt`, `omega`, `vref`,
     * `vs`, `vuel` and `voel` of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Ieeet1<scalar_type, index_type>::Ieeet1(const ModelDataT& data)
      : ComponentT(std::type_identity<Ieeet1>{}),
        V_(data.V),
        Tr_(std::max(data.Tr, TIME_CONSTANT_MINIMUM)),
        Ka_(data.Ka),
        Ta_(std::max(data.Ta, TIME_CONSTANT_MINIMUM)),
        Ke_(data.Ke),
        Te_(std::max(data.Te, TIME_CONSTANT_MINIMUM)),
        Kf_(data.Kf),
        Tf_(std::max(data.Tf, TIME_CONSTANT_MINIMUM)),
        Vrmin_(data.Vrmin),
        Vrmax_(data.Vrmax),
        E1_(data.E1),
        E2_(data.E2),
        Se1_(data.Se1),
        Se2_(data.Se2),
        Ispdlim_(data.Ispdlim)
    {
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
     * speed multiplier flag and the saturation data.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Ieeet1<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Ieeet1");
      checks.check(V_ > ZERO<RealT>, "V must be positive");
      checks.check(Ka_ > ZERO<RealT>, "Ka must be positive");
      checks.check(Ke_ != ZERO<RealT>, "Ke must be nonzero: a value from the operating point is not supported");
      checks.check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
      checks.check(Ispdlim_ == ZERO<RealT> || Ispdlim_ == ONE<RealT>, "Ispdlim must be 0 or 1");
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
     * @brief Write the residuals of the voltage sensing, of the regulator, of
     * the exciter and of the stabilizing feedback.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Ieeet1<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T vt    = w[VT] / V_;
      const T omega = w[OMEGA] - ONE<RealT>;

      f[VTS]  = -yp[VTS] + (vt - y[VTS]) / Tr_;
      f[VR]   = -yp[VR] + Math::antiwindup(y[VR], (-y[VR] + Ka_ * y[VTR]) / Ta_, Vrmin_, Vrmax_);
      f[EFDP] = -yp[EFDP] + (y[VR] - y[VE] - Ke_ * y[EFDP]) / Te_;
      f[VFX]  = -yp[VFX] + y[VF] / Tf_;
      f[VTR]  = -y[VTS] + w[VREF] + w[VS] + w[VUEL] + w[VOEL] - y[VTR] - y[VF];
      f[VF]   = -Tf_ * (y[VF] + y[VFX]) + Kf_ * y[EFDP];
      f[VE]   = -y[VE] + y[KSAT];
      f[EFD]  = -y[EFD] + y[EFDP] + omega * y[EFDP] * Ispdlim_;
      f[KSAT] = -y[KSAT] + SB_ * Math::qramp(y[EFDP] - SA_);
    }
  } // namespace EMT
} // namespace GridKit
