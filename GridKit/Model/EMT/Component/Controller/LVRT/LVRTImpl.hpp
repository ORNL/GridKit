/**
 * @file LVRTImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of EMT low-voltage ride-through.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/LVRT/LVRT.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a low-voltage ride-through.
     *
     * Layout: internal variables `il` and `ip` of size 1, and inputs `vt`
     * and `ipcmd` of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    LVRT<scalar_type, index_type>::LVRT(const ModelDataT& data)
      : ComponentT(std::type_identity<LVRT>{}),
        TM_(data.TM),
        VL0_(data.VL0),
        VL1_(data.VL1),
        IL1_(data.IL1),
        Rpmax_(data.Rpmax),
        VA0_(data.VA0),
        VA1_(data.VA1)
    {
    }

    /**
     * @brief Check the time constant, the voltage breakpoints, the ceiling
     *        target and the recovery rate.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks LVRT<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("LVRT");
      checks.check(TM_ > ZERO<RealT>, "TM must be positive");
      checks.check(VL0_ >= ZERO<RealT> && VL1_ > VL0_, "VL0 must be nonnegative and below VL1");
      checks.check(VA0_ >= ZERO<RealT> && VA1_ > VA0_, "VA0 must be nonnegative and below VA1");
      checks.check(IL1_ > ZERO<RealT> && Rpmax_ > ZERO<RealT>, "IL1 and Rpmax must be positive");
      return checks;
    }

    /**
     * @brief Write the residuals of the current ceiling and of the limited
     *        active-current order.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    LVRT<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto IL     = static_cast<size_t>(LVRTInternalVariables::il);
      const auto IP     = static_cast<size_t>(LVRTInternalVariables::ip);
      const auto VT     = static_cast<size_t>(LVRTInputs::vt);
      const auto IP_CMD = static_cast<size_t>(LVRTInputs::ipcmd);

      const T target    = Math::linseg(w[VT], VL0_, VL1_, IL1_);
      const T rate      = (target - y[IL]) / TM_;
      const T available = y[IL] * Math::linseg(w[VT], VA0_, VA1_, ONE<RealT>);

      f[IL] = -yp[IL] + std::min(rate, T{Rpmax_});
      f[IP] = -y[IP] + w[IP_CMD] - std::fdim(w[IP_CMD], available);
    }
  } // namespace EMT
} // namespace GridKit
