/**
 * @file Esst1aImpl.hpp
 * @brief Equations of the IEEE ST1A static excitation system.
 */

#pragma once

#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esst1a/Esst1a.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    Esst1a<ScalarT, IdxT>::Esst1a(const ModelDataT& data)
      : ComponentT(std::type_identity<Esst1a>{}), data_(data)
    {
      washout_gain_     = data.Tf > 0 ? data.Kf / data.Tf : 0;
      lead_gain_        = data.Tb > 0 ? data.Tc / data.Tb : 0;
      second_lead_gain_ = data.Tb1 > 0 ? data.Tc1 / data.Tb1 : 0;
    }

    template <typename ScalarT, typename IdxT>
    Model::ConfigurationChecks Esst1a<ScalarT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Esst1a");
      const auto&                d = data_;
      checks.check(d.V > 0 && d.Ka > 0, "V and Ka must be positive");
      checks.check(d.Tr >= 0 && d.Ta >= 0, "Tr and Ta must be nonnegative");
      checks.check(d.Tc >= 0 && d.Tb >= 0 && (d.Tb > 0 || d.Tc == 0), "Tc requires a positive Tb");
      checks.check(d.Tf >= 0 && (d.Tf > 0 || d.Kf == 0), "Kf requires a positive Tf");
      checks.check(d.Vimin <= d.Vimax && d.Vrmin <= d.Vrmax, "voltage limits must be ordered");
      checks.check(d.Tc1 >= 0 && d.Tb1 >= 0 && (d.Tb1 > 0 || d.Tc1 == 0), "Tc1 requires a positive Tb1");
      checks.check(d.Vamin <= d.Vamax, "regulator limits must be ordered");
      checks.check(d.UEL >= 0 && d.UEL <= 3 && (d.VOS == 1 || d.VOS == 2), "invalid UEL or VOS input location");
      return checks;
    }

    template <typename ScalarT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] void
    Esst1a<ScalarT, IdxT>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto& d       = data_;
      const T     vt      = w[VT_IN] / d.V;
      const T     ifd     = w[IFD_IN];
      const T     limited = y[VA] + (d.VOS == 2) * w[VS_IN] - d.Klr * Math::ramp(ifd - d.Ilr);
      const T     gated   = limited + (d.UEL == 3) * Math::ramp(w[VUEL_IN] - limited);
      f[REG]              = Math::min(gated, w[VOEL_IN]) - y[REG];
      const T feedback    = washout_gain_ * (y[REG] - y[WASH]);
      const T error       = w[VREF_IN] - y[VC] - feedback + (d.VOS == 1) * w[VS_IN]
                      + (d.UEL == 1) * w[VUEL_IN];
      const T vi = Math::clamp(error, d.Vimin, d.Vimax);
      f[DRIVE]   = vi + (d.UEL == 2) * Math::ramp(w[VUEL_IN] - vi) - y[DRIVE];
      f[STAGE]   = y[LL] + lead_gain_ * (y[DRIVE] - y[LL]) - y[STAGE];
      f[DRIVER]  = d.Ka * (y[LL1] + second_lead_gain_ * (y[STAGE] - y[LL1])) - y[DRIVER];

      f[VC]  = vt - y[VC] - d.Tr * yp[VC];
      f[LL]  = y[DRIVE] - y[LL] - d.Tb * yp[LL];
      f[LL1] = y[STAGE] - y[LL1] - d.Tb1 * yp[LL1];
      f[VA]  = (d.Ta > 0) * Math::antiwindup(y[VA], y[DRIVER] - y[VA], d.Vamin, d.Vamax)
              + (d.Ta == 0) * (Math::clamp(y[DRIVER], d.Vamin, d.Vamax) - y[VA]) - d.Ta * yp[VA];
      f[WASH] = y[REG] - y[WASH] - d.Tf * yp[WASH];
      f[EFD]  = Math::max(vt * d.Vrmin, Math::min(y[REG], vt * d.Vrmax - d.Kc * ifd)) - y[EFD];
    }
  } // namespace EMT
} // namespace GridKit
