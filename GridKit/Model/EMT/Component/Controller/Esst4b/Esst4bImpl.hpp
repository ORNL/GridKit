/**
 * @file Esst4bImpl.hpp
 * @brief Equations of the IEEE ST4B static excitation system.
 */

#pragma once

#include <cmath>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Esst4b/Esst4b.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/ThreePhase.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename ScalarT, typename IdxT>
    Esst4b<ScalarT, IdxT>::Esst4b(const ModelDataT& data)
      : ComponentT(std::type_identity<Esst4b>{}), data_(data)
    {
      if (data.V > 0 && data.S > 0)
      {
        inv_v_peak_ = std::sqrt(1.5) / data.V;
        inv_i_peak_ = std::sqrt(1.5) * data.V / data.S;
      }
      const RealT angle = data.ThetaP * std::numbers::pi / 180;
      kp_real_          = data.Kp * std::cos(angle);
      kp_imag_          = data.Kp * std::sin(angle);
    }

    template <typename ScalarT, typename IdxT>
    Model::ConfigurationChecks Esst4b<ScalarT, IdxT>::verify() const
    {
      Model::ConfigurationChecks checks("Esst4b");
      const auto&                d = data_;
      checks.check(d.V > 0 && d.S > 0, "V and S must be positive");
      checks.check(d.Tr >= 0 && d.Ta >= 0, "Tr and Ta must be nonnegative");
      checks.check(d.Vrmin <= d.Vrmax && d.Vmmin <= d.Vmmax, "regulator limits must be ordered");
      checks.check(d.Kc >= 0 && d.Vbmax > 0, "Kc must be nonnegative and Vbmax positive");
      return checks;
    }

    template <typename ScalarT, typename IdxT>
    template <typename T>
    [[gnu::always_inline]] void
    Esst4b<ScalarT, IdxT>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto& d       = data_;
      // Clarke components in peak per unit. The measured current is positive into the grid.
      const auto  voltage = clarke(w + V_IN);
      const auto  current = clarke(w + I_IN);
      const T     vr      = voltage.alpha * inv_v_peak_;
      const T     vi      = voltage.beta * inv_v_peak_;
      const T     ir      = current.alpha * inv_i_peak_;
      const T     ii      = current.beta * inv_i_peak_;
      const RealT a = kp_real_, b = kp_imag_;
      const T     er     = a * vr - b * vi - d.Xl * b * ir - (d.Ki + d.Xl * a) * ii;
      const T     ei     = b * vr + a * vi + (d.Ki + d.Xl * a) * ir - d.Xl * b * ii;
      const T     square = er * er + ei * ei;
      // V_E F_EX, evaluated on its three loading ranges without dividing by V_E.
      const T     ve     = square / std::sqrt(std::max(square, T{1e-24}));
      const T     load   = d.Kc * w[IFD_IN];
      const T     knee   = (std::numbers::sqrt3 / 4) * ve;
      const T     first  = std::fdim(load, T{0}) - std::fdim(load, knee);
      const T     middle = knee + std::fdim(load, knee) - std::fdim(load, 0.75 * ve);
      const T     last   = 0.75 * ve + std::fdim(load, 0.75 * ve) - std::fdim(load, ve);
      const T     arc    = 0.75 * square - middle * middle;
      const T     bridge = 0.25 * ve - first / std::numbers::sqrt3 + arc / std::sqrt(std::max(arc, T{1e-24}))
                       - std::numbers::sqrt3 * (last - 0.75 * ve);

      const T error        = w[VREF_IN] + w[VS_IN] + w[VUEL_IN] - y[VC];
      const T vr_unlimited = d.Kpr * error + y[XVR];
      const T field_error  = y[VA] - Math::min(d.Kg * y[EFD], d.Vgmax);
      const T vm           = d.Kpm * field_error + y[XVM];

      f[VC]  = w[VT_IN] / d.V - y[VC] - d.Tr * yp[VC];
      f[XVR] = Math::antiwindup(vr_unlimited, d.Kir * error, d.Vrmin, d.Vrmax) - yp[XVR];
      f[VA]  = Math::clamp(vr_unlimited, d.Vrmin, d.Vrmax) - y[VA] - d.Ta * yp[VA];
      f[XVM] = Math::antiwindup(vm, d.Kim * field_error, d.Vmmin, d.Vmmax) - yp[XVM];
      f[EFD] = Math::min(bridge, d.Vbmax) * Math::min(Math::clamp(vm, d.Vmmin, d.Vmmax), w[VOEL_IN]) - y[EFD];
    }
  } // namespace EMT
} // namespace GridKit
