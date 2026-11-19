/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT IEEE power system stabilizer.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeest/Ieeest.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct an IEEE power system stabilizer and the selections of its blocks.
     *
     * Layout: 12 internal variables of size 1, and input `u` of size 1.
     * `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Ieeest<scalar_type, index_type>::Ieeest(const ModelDataT& data)
      : ComponentT(std::type_identity<Ieeest>{}),
        A5_(data.A5),
        A6_(data.A6),
        T1_(data.T1),
        T2_(data.T2),
        T3_(data.T3),
        T4_(data.T4),
        T5_(data.T5),
        T6_(data.T6),
        Ks_(data.Ks),
        Lsmin_(data.Lsmin),
        Lsmax_(data.Lsmax),
        Vcl_(data.Vcl),
        Vcu_(data.Vcu),
        Tdelay_(data.Tdelay)
    {
      a1_ = data.A1 + data.A3;
      a2_ = data.A2 + data.A4 + data.A1 * data.A3;
      a3_ = data.A1 * data.A4 + data.A2 * data.A3;
      a4_ = data.A2 * data.A4;

      // The residual has no branches, so the order of the notch filter and each bypass are selections
      use_4th_order_ = static_cast<RealT>(a4_ != ZERO<RealT>);
      use_3rd_order_ = static_cast<RealT>(a4_ == ZERO<RealT> && a3_ != ZERO<RealT>);
      use_2nd_order_ = static_cast<RealT>(a4_ == ZERO<RealT> && a3_ == ZERO<RealT> && a2_ != ZERO<RealT>);
      use_notch_     = use_4th_order_ + use_3rd_order_ + use_2nd_order_;
      inv_a4_        = use_4th_order_ / (a4_ + (ONE<RealT> - use_4th_order_));
      inv_a3_        = use_3rd_order_ / (a3_ + (ONE<RealT> - use_3rd_order_));
      inv_a2_        = use_2nd_order_ / (a2_ + (ONE<RealT> - use_2nd_order_));
      use_T2_        = static_cast<RealT>(T2_ != ZERO<RealT>);
      use_T4_        = static_cast<RealT>(T4_ != ZERO<RealT>);
      use_T6_        = static_cast<RealT>(T6_ != ZERO<RealT>);
    }

    /**
     * @brief Check the order of the notch filter and the parameters that are
     * not supported.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Ieeest<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Ieeest");
      checks.check(use_notch_ == ONE<RealT> || a1_ == ZERO<RealT>, "a first-order notch filter is not supported");
      checks.check(Lsmin_ <= Lsmax_, "Lsmin must be less than or equal to Lsmax");
      checks.check(Vcl_ == ZERO<RealT> && Vcu_ == ZERO<RealT>, "Vcl and Vcu must be zero: not supported");
      checks.check(Tdelay_ == ZERO<RealT>, "Tdelay must be zero: not supported");
      return checks;
    }

    /**
     * @brief Write the residuals of the notch filter, of the lead-lag blocks,
     * of the washout and of the limited signal.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Ieeest<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T u      = w[U];
      const T x2_rhs = (use_4th_order_ + use_3rd_order_) * y[X3]
                       + use_2nd_order_ * (-y[X1] - a1_ * y[X2] + u) * inv_a2_;

      f[X1] = -yp[X1] + use_notch_ * y[X2];
      f[X2] = -yp[X2] + x2_rhs;
      f[X3] = -yp[X3] + use_4th_order_ * y[X4]
              + use_3rd_order_ * (-y[X1] - a1_ * y[X2] - a2_ * y[X3] + u) * inv_a3_;
      f[X4]  = -yp[X4] + use_4th_order_ * (-y[X1] - a1_ * y[X2] - a2_ * y[X3] - a3_ * y[X4] + u) * inv_a4_;
      f[X5]  = -T2_ * yp[X5] - y[X5] + y[V4];
      f[X6]  = -T4_ * yp[X6] - y[X6] + y[V5];
      f[X7]  = -T6_ * yp[X7] - y[X7] + y[V6];
      f[V4]  = -y[V4] + (ONE<RealT> - use_notch_) * u + use_notch_ * (y[X1] + A5_ * y[X2] + A6_ * x2_rhs);
      f[V5]  = use_T2_ * (-T2_ * (y[V5] - y[X5]) + T1_ * (y[V4] - y[X5])) + (ONE<RealT> - use_T2_) * (y[V4] - y[V5]);
      f[V6]  = use_T4_ * (-T4_ * (y[V6] - y[X6]) + T3_ * (y[V5] - y[X6])) + (ONE<RealT> - use_T4_) * (y[V5] - y[V6]);
      f[V7]  = use_T6_ * (-T6_ * y[V7] + Ks_ * T5_ * (y[V6] - y[X7])) + (ONE<RealT> - use_T6_) * (Ks_ * y[V6] - y[V7]);
      f[VSS] = -y[VSS] + Math::clamp(y[V7], Lsmin_, Lsmax_);
    }
  } // namespace EMT
} // namespace GridKit
