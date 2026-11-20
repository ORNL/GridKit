/**
 * @file SexsPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT simplified excitation system.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/SexsPti/SexsPti.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a simplified excitation system.
     *
     * Layout: 3 internal variables of size 1, and inputs `vt`, `vref`, `vs`,
     * `vuel` and `voel` of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    SexsPti<scalar_type, index_type>::SexsPti(const ModelDataT& data)
      : ComponentT(std::type_identity<SexsPti>{}),
        Ta_(data.Ta),
        Tb_(data.Tb),
        Te_(data.Te),
        K_(data.K),
        Efdmax_(data.Efdmax),
        Efdmin_(data.Efdmin),
        V_(data.V)
    {
    }

    /**
     * @brief Check the rated voltage, the time constants, the gain and the
     * field voltage limits.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks SexsPti<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("SexsPti");
      checks.check(V_ > ZERO<RealT>, "V must be positive");
      checks.check(Ta_ >= ZERO<RealT>, "Ta must be nonnegative");
      checks.check(Tb_ > ZERO<RealT>, "Tb must be positive");
      checks.check(Te_ > ZERO<RealT>, "Te must be positive");
      checks.check(K_ > ZERO<RealT>, "K must be positive");
      checks.check(Efdmin_ < Efdmax_, "Efdmin must be less than Efdmax");
      return checks;
    }

    /**
     * @brief Write the residuals of the lead-lag block, of the field voltage
     * and of the voltage error.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    SexsPti<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T vt   = w[VT] / V_;
      const T lead = -y[VR] + Ta_ * y[VTR];

      f[VR]  = -yp[VR] + lead / Tb_ - y[VTR];
      f[EFD] = -yp[EFD] + Math::antiwindup(y[EFD], (-y[EFD] + (K_ / Tb_) * lead) / Te_, Efdmin_, Efdmax_);
      f[VTR] = -y[VTR] - vt + w[VREF] + w[VS] + w[VUEL] + w[VOEL];
    }
  } // namespace EMT
} // namespace GridKit
