/**
 * @file GastPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT gas turbine governor.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/GastPti/GastPti.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a gas turbine governor and its rating ratio.
     *
     * Layout: 7 internal variables of size 1, and inputs `omega` and `pref` of
     * size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    GastPti<scalar_type, index_type>::GastPti(const ModelDataT& data)
      : ComponentT(std::type_identity<GastPti>{}),
        R_(data.R),
        T1_(std::max(data.T1, TIME_CONSTANT_MINIMUM)),
        T2_(std::max(data.T2, TIME_CONSTANT_MINIMUM)),
        T3_(std::max(data.T3, TIME_CONSTANT_MINIMUM)),
        At_(data.At),
        Kt_(data.Kt),
        Vmax_(data.Vmax),
        Vmin_(data.Vmin),
        Dturb_(data.Dturb),
        Trate_(data.Trate),
        S_(data.S)
    {
      if (Trate_ > ZERO<RealT>)
      {
        machine_to_turbine_ = S_ / (MEGA<RealT> * Trate_);
      }
      if (Vmin_ < Vmax_)
      {
        valve_on_ = ONE<RealT>;
      }
    }

    /**
     * @brief Check the ratings, the droop, the temperature limit and its gain,
     * the valve limits and the damping.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks GastPti<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("GastPti");
      checks.check(Trate_ > ZERO<RealT>, "Trate must be positive");
      checks.check(S_ > ZERO<RealT>, "S must be positive");
      checks.check(std::isfinite(R_) && R_ > ZERO<RealT>, "R must be finite and positive");
      checks.check(std::isfinite(At_) && At_ >= ZERO<RealT>, "At must be finite and nonnegative");
      checks.check(std::isfinite(Kt_) && Kt_ >= ZERO<RealT>, "Kt must be finite and nonnegative");
      checks.check(std::isfinite(Vmin_) && std::isfinite(Vmax_) && Vmin_ <= Vmax_,
                   "Vmin and Vmax must be finite, with Vmin less than or equal to Vmax");
      checks.check(std::isfinite(Dturb_) && Dturb_ >= ZERO<RealT>, "Dturb must be finite and nonnegative");
      return checks;
    }

    /**
     * @brief Write the residuals of the fuel system, of the fuel demands and of
     * the mechanical power.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    GastPti<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const RealT system_to_turbine = this->va_system_base_ / (MEGA<RealT> * Trate_);

      const T pmech = machine_to_turbine_ * y[PMECH];
      const T omega = w[OMEGA] - ONE<RealT>;
      const T pref  = system_to_turbine * w[PREF];

      f[XVALVE] = -yp[XVALVE] + valve_on_ * Math::antiwindup(y[XVALVE], y[VLV] - y[XVALVE], Vmin_, Vmax_) / T1_;
      f[XFLOW]  = -yp[XFLOW] + (-y[XFLOW] + y[XVALVE]) / T2_;
      f[XTEMP]  = -yp[XTEMP] + (-y[XTEMP] + y[XFLOW]) / T3_;
      f[VLOAD]  = -omega + R_ * (pref - y[VLOAD]);
      f[VTEMP]  = -y[VTEMP] + At_ + Kt_ * (At_ - y[XTEMP]);
      f[VLV]    = -y[VLV] + Math::min(y[VLOAD], y[VTEMP]);
      f[PMECH]  = -pmech + y[XFLOW] - Dturb_ * omega;
    }
  } // namespace EMT
} // namespace GridKit
