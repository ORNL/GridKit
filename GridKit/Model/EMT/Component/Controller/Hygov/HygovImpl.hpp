/**
 * @file HygovImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT hydro turbine governor.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Hygov/Hygov.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a hydro turbine governor, its gain curve and its rating ratio.
     *
     * Layout: 12 internal variables of size 1, and inputs `omega`, `pref` and
     * `paux` of size 1. A gain curve of zeros selects the default curve.
     * `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Hygov<scalar_type, index_type>::Hygov(const ModelDataT& data)
      : ComponentT(std::type_identity<Hygov>{}),
        Rperm_(data.Rperm),
        Rtemp_(data.Rtemp),
        Tr_(std::max(data.Tr, TIME_CONSTANT_MINIMUM)),
        Tf_(std::max(data.Tf, TIME_CONSTANT_MINIMUM)),
        Tg_(std::max(data.Tg, TIME_CONSTANT_MINIMUM)),
        Velm_(data.Velm),
        Gmax_(data.Gmax),
        Gmin_(data.Gmin),
        Tw_(std::max(data.Tw, TIME_CONSTANT_MINIMUM)),
        At_(data.At),
        Dturb_(data.Dturb),
        Qnl_(data.Qnl),
        Tn_(data.Tn),
        Tnp_(std::max(data.Tnp, TIME_CONSTANT_MINIMUM)),
        db1_(data.db1),
        db2_(data.db2),
        Hdam_(data.Hdam),
        Trate_(data.Trate),
        S_(data.S)
    {
      Gv_  = {data.Gv0, data.Gv1, data.Gv2, data.Gv3, data.Gv4, data.Gv5};
      Pgv_ = {data.Pgv0, data.Pgv1, data.Pgv2, data.Pgv3, data.Pgv4, data.Pgv5};

      const auto nonzero = [](RealT value)
      { return value != ZERO<RealT>; };

      if (std::ranges::none_of(Gv_, nonzero) && std::ranges::none_of(Pgv_, nonzero))
      {
        Gv_  = DEFAULT_CURVE;
        Pgv_ = DEFAULT_CURVE;
      }

      leadlag_gain_ = Tn_ / Tnp_;
      if (Trate_ > ZERO<RealT>)
      {
        machine_to_turbine_ = S_ / (MEGA<RealT> * Trate_);
      }
    }

    /**
     * @brief Check the ratings, the droops, the limits, the deadbands, the head
     * and the gain curve.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Hygov<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Hygov");
      checks.check(Trate_ > ZERO<RealT>, "Trate must be positive");
      checks.check(S_ > ZERO<RealT>, "S must be positive");
      checks.check(Rtemp_ > ZERO<RealT>, "Rtemp must be positive");
      checks.check(Tn_ >= ZERO<RealT>, "Tn must be nonnegative");
      checks.check(Velm_ >= ZERO<RealT>, "Velm must be nonnegative");
      checks.check(Gmin_ < Gmax_, "Gmin must be less than Gmax");
      checks.check(At_ > ZERO<RealT>, "At must be positive");
      checks.check(Dturb_ >= ZERO<RealT>, "Dturb must be nonnegative");
      checks.check(db1_ >= ZERO<RealT>, "db1 must be nonnegative");
      checks.check(db2_ == ZERO<RealT>, "db2 must be zero: backlash is not supported");
      checks.check(Hdam_ > ZERO<RealT>, "Hdam must be positive");
      for (size_t k = 1; k < CURVE_POINTS; ++k)
      {
        checks.check(Gv_[k - 1] < Gv_[k], "Gv points must be strictly increasing");
        checks.check(Pgv_[k - 1] <= Pgv_[k], "Pgv points must be nondecreasing");
      }
      checks.check(Gv_.front() <= Gmin_, "Gmin must be at or above the first Gv point");
      checks.check(Gmax_ <= Gv_.back(), "Gmax must be at or below the last Gv point");
      return checks;
    }

    /**
     * @brief Power of the gain curve at a gate position.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] gate Gate position.
     * @return The power.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Hygov<scalar_type, index_type>::gatePower(T gate) const
    {
      return Pgv_[0]
             + Math::linseg(gate, Gv_[0], Gv_[1], Pgv_[1] - Pgv_[0])
             + Math::linseg(gate, Gv_[1], Gv_[2], Pgv_[2] - Pgv_[1])
             + Math::linseg(gate, Gv_[2], Gv_[3], Pgv_[3] - Pgv_[2])
             + Math::linseg(gate, Gv_[3], Gv_[4], Pgv_[4] - Pgv_[3])
             + Math::linseg(gate, Gv_[4], Gv_[5], Pgv_[5] - Pgv_[4]);
    }

    /**
     * @brief Write the residuals of the governor, of the gate, of the penstock
     * and of the mechanical power.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Hygov<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const RealT system_to_turbine = this->va_system_base_ / (MEGA<RealT> * Trate_);

      const T pmech  = machine_to_turbine_ * y[PMECH];
      const T omega  = w[OMEGA] - ONE<RealT>;
      const T pref   = system_to_turbine * (w[PREF] + w[PAUX]);
      const T yomega = y[XN] + leadlag_gain_ * (y[OMEGADB] - y[XN]);

      f[XN]      = -yp[XN] + (y[OMEGADB] - y[XN]) / Tnp_;
      f[XF]      = -yp[XF] + (y[EF] - y[XF]) / Tf_;
      f[C]       = -yp[C] + Math::antiwindup(y[C], y[RC], Gmin_, Gmax_);
      f[G]       = -yp[G] + (y[C] - y[G]) / Tg_;
      f[Q]       = -yp[Q] + (Hdam_ - y[H]) / Tw_;
      f[OMEGADB] = -y[OMEGADB] + Math::deadband1(omega, -db1_, db1_);
      f[EF]      = -y[EF] + pref - yomega - Rperm_ * y[C];
      f[FC]      = -Rtemp_ * y[FC] + y[XF] / Tr_ + (y[EF] - y[XF]) / Tf_;
      f[RC]      = -y[RC] + Math::clamp(y[FC], -Velm_, Velm_);
      f[PGV]     = -y[PGV] + gatePower(y[G]);
      f[H]       = -y[Q] * y[Q] + y[H] * y[PGV] * y[PGV];
      f[PMECH]   = -pmech + At_ * y[H] * (y[Q] - Qnl_) - Dturb_ * omega * y[G];
    }
  } // namespace EMT
} // namespace GridKit
