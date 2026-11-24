/**
 * @file Exst1PtiImpl.hpp
 * @brief Equations of the PTI static excitation system.
 */

#pragma once

#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Exst1Pti/Exst1Pti.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    Exst1Pti<ScalarT, IdxT>::Exst1Pti(const ModelDataT& data)
      : ComponentT(std::type_identity<Exst1Pti>{}), data_(data)
    {
      washout_gain_ = data.Tf > 0 ? data.Kf / data.Tf : 0;
      lead_gain_    = data.Tb > 0 ? data.Tc / data.Tb : 0;
    }

    template <typename ScalarT, typename IdxT>
    Model::ConfigurationChecks Exst1Pti<ScalarT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Exst1Pti");
      const auto&                d = data_;
      checks.check(d.V > 0 && d.Ka > 0, "V and Ka must be positive");
      checks.check(d.Tr >= 0 && d.Ta >= 0, "Tr and Ta must be nonnegative");
      checks.check(d.Tc >= 0 && d.Tb >= 0 && (d.Tb > 0 || d.Tc == 0), "Tc requires a positive Tb");
      checks.check(d.Tf >= 0 && (d.Tf > 0 || d.Kf == 0), "Kf requires a positive Tf");
      checks.check(d.Vimin <= d.Vimax && d.Vrmin <= d.Vrmax, "voltage limits must be ordered");
      return checks;
    }

    template <typename ScalarT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] void
    Exst1Pti<ScalarT, IdxT>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto& d        = data_;
      const T     vt       = w[VT_IN] / d.V;
      const T     feedback = washout_gain_ * (y[VA] - y[WASH]);
      const T     error    = Math::clamp(w[VREF_IN] + w[VS_IN] + w[VUEL_IN] + w[VOEL_IN] - y[VC] - feedback, d.Vimin, d.Vimax);
      const T     lead     = y[LL] + lead_gain_ * (error - y[LL]);
      const T     drop     = d.Kc * w[IFD_IN];

      f[VC]   = vt - y[VC] - d.Tr * yp[VC];
      f[LL]   = error - y[LL] - d.Tb * yp[LL];
      f[VA]   = d.Ka * lead - y[VA] - d.Ta * yp[VA];
      f[WASH] = y[VA] - y[WASH] - d.Tf * yp[WASH];
      f[EFD]  = Math::clamp(y[VA], vt * d.Vrmin - drop, vt * d.Vrmax - drop) - y[EFD];
    }
  } // namespace EMT
} // namespace GridKit
