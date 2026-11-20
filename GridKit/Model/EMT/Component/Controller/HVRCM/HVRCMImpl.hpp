/**
 * @file HVRCMImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of EMT high-voltage reactive current management.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/HVRCM/HVRCM.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a high-voltage reactive current management.
     *
     * Layout: internal variable `iq` of size 1, and inputs `vt` and `iqcmd`
     * of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    HVRCM<scalar_type, index_type>::HVRCM(const ModelDataT& data)
      : ComponentT(std::type_identity<HVRCM>{}),
        Vhvmax_(data.Vhvmax),
        Khv_(data.Khv),
        Qmin_(data.Qmin)
    {
    }

    /**
     * @brief Check the voltage threshold, the gain and the lower bound.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks HVRCM<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("HVRCM");
      checks.check(Vhvmax_ > ZERO<RealT>, "Vhvmax must be positive");
      checks.check(Khv_ >= ZERO<RealT>, "Khv must be nonnegative");
      checks.check(Qmin_ <= ZERO<RealT>, "Qmin must be nonpositive");
      return checks;
    }

    /**
     * @brief Write the residual of the reduced reactive-current order.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    HVRCM<scalar_type, index_type>::internalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const auto IQ     = static_cast<size_t>(HVRCMInternalVariables::iq);
      const auto VT     = static_cast<size_t>(HVRCMInputs::vt);
      const auto IQ_CMD = static_cast<size_t>(HVRCMInputs::iqcmd);

      const T absorbed = Khv_ * Math::ramp(w[VT] - Vhvmax_);
      const T floor    = std::min(w[IQ_CMD], T{Qmin_});

      f[IQ] = -y[IQ] + floor + std::fdim(w[IQ_CMD] - floor, absorbed);
    }
  } // namespace EMT
} // namespace GridKit
