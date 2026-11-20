/**
 * @file Tgov1Impl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT steam turbine governor.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Tgov1/Tgov1.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a steam turbine governor and its rating ratio.
     *
     * Layout: 3 internal variables of size 1, and inputs `omega` and `pref` of
     * size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Tgov1<scalar_type, index_type>::Tgov1(const ModelDataT& data)
      : ComponentT(std::type_identity<Tgov1>{}),
        Trate_(data.Trate),
        R_(data.R),
        Pvmin_(data.Pvmin),
        Pvmax_(data.Pvmax),
        T1_(std::max(data.T1, TIME_CONSTANT_MINIMUM)),
        T2_(data.T2),
        T3_(std::max(data.T3, TIME_CONSTANT_MINIMUM)),
        Dt_(data.Dt),
        S_(data.S)
    {
      if (Trate_ > ZERO<RealT>)
      {
        machine_to_turbine_ = S_ / (MEGA<RealT> * Trate_);
      }
    }

    /**
     * @brief Check the ratings, the droop and the valve limits.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Tgov1<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Tgov1");
      checks.check(Trate_ > ZERO<RealT>, "Trate must be positive");
      checks.check(S_ > ZERO<RealT>, "S must be positive");
      checks.check(R_ != ZERO<RealT>, "R must be nonzero");
      checks.check(Pvmin_ <= Pvmax_, "Pvmin must be less than or equal to Pvmax");
      return checks;
    }

    /**
     * @brief Write the residuals of the turbine, of the valve and of the
     * mechanical power.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Tgov1<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T pm    = machine_to_turbine_ * y[PM];
      const T omega = w[OMEGA] - ONE<RealT>;

      f[PTX] = -yp[PTX] - (y[PTX] - y[PV] - T2_ * yp[PV]) / T3_;
      f[PV]  = -yp[PV] + Math::antiwindup(y[PV], -y[PV] + (w[PREF] - omega) / R_, Pvmin_, Pvmax_) / T1_;
      f[PM]  = -pm + y[PTX] - Dt_ * omega;
    }
  } // namespace EMT
} // namespace GridKit
