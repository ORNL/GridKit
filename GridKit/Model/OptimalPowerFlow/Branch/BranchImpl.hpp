/**
 * @file BranchImpl.hpp
 * @brief Definition of an optimal power flow branch.
 */

#pragma once

#include <array>
#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::Branch(const ModelDataT& data)
      : ComponentT(std::type_identity<Branch>{}, data.id),
        in_service_(data.in_service),
        R_(parameter(data, BranchParameters::R, ZERO<RealT>)),
        X_(parameter(data, BranchParameters::X, ZERO<RealT>)),
        G_(parameter(data, BranchParameters::G, ZERO<RealT>)),
        B_(parameter(data, BranchParameters::B, ZERO<RealT>)),
        Gmag_(parameter(data, BranchParameters::Gmag, ZERO<RealT>)),
        Bmag_(parameter(data, BranchParameters::Bmag, ZERO<RealT>)),
        tap_(parameter(data, BranchParameters::tap, ONE<RealT>)),
        phase_(parameter(data, BranchParameters::phase, ZERO<RealT>)),
        Smax_(parameter(data, BranchParameters::Smax, UNBOUNDED<RealT>))
    {
      setDerivedParams();
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Branch<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Branch");
      checks.check(std::isfinite(R_), "R must be finite");
      checks.check(std::isfinite(X_), "X must be finite");
      checks.check(std::isfinite(G_), "G must be finite");
      checks.check(std::isfinite(B_), "B must be finite");
      checks.check(std::isfinite(Gmag_), "Gmag must be finite");
      checks.check(std::isfinite(Bmag_), "Bmag must be finite");
      checks.check(std::isfinite(tap_), "tap must be finite");
      checks.check(std::isfinite(phase_), "phase must be finite");
      checks.check(R_ * R_ + X_ * X_ > ZERO<RealT>, "R and X cannot both be zero");
      checks.check(tap_ > ZERO<RealT>, "tap must be positive");
      checks.check(Smax_ > ZERO<RealT>, "Smax must be positive");
      return checks;
    }

    /**
     * @brief Apparent power at most the rating at each port
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::allocate()
    {
      RealT* upper = d_upper_.getData();
      upper[0]     = Smax_ * Smax_;
      upper[1]     = Smax_ * Smax_;
      return 0;
    }

    /**
     * @brief Active and reactive power into bus 1, then into bus 2, added to
     * their balances
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Branch<scalar_type, index_type>::externalEqualities(const ScalarT* x, ScalarT* c) const
    {
      power(x, c);
    }

    /**
     * @brief Squared apparent power at bus 1 and at bus 2, limited by the
     * rating
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Branch<scalar_type, index_type>::internalInequalities(const ScalarT* x, ScalarT* d) const
    {
      std::array<ScalarT, 4> s;
      power(x, s.data());

      d[0] = s[0] * s[0] + s[1] * s[1];
      d[1] = s[2] * s[2] + s[3] * s[3];
    }

    /**
     * @brief Active and reactive power into bus 1, then into bus 2
     *
     * The power into each bus is \f$P + jQ = V I^*\f$.
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Branch<scalar_type, index_type>::power(const ScalarT* x, ScalarT* s) const
    {
      const ScalarT vr1 = x[0];
      const ScalarT vi1 = x[1];
      const ScalarT vr2 = x[2];
      const ScalarT vi2 = x[3];

      const ScalarT ir1 = g11_ * vr1 - b11_ * vi1 + g12_ * vr2 - b12_ * vi2;
      const ScalarT ii1 = b11_ * vr1 + g11_ * vi1 + b12_ * vr2 + g12_ * vi2;
      const ScalarT ir2 = g21_ * vr1 - b21_ * vi1 + g22_ * vr2 - b22_ * vi2;
      const ScalarT ii2 = b21_ * vr1 + g21_ * vi1 + b22_ * vr2 + g22_ * vi2;

      s[0] = vr1 * ir1 + vi1 * ii1;
      s[1] = vi1 * ir1 - vr1 * ii1;
      s[2] = vr2 * ir2 + vi2 * ii2;
      s[3] = vi2 * ir2 - vr2 * ii2;
    }

    /**
     * @brief Pi-model admittance with the off-nominal tap and phase shift on
     * the bus-1 side
     */
    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::setDerivedParams()
    {
      const RealT denom   = R_ * R_ + X_ * X_;
      const RealT g_br    = R_ / denom;
      const RealT b_br    = -X_ / denom;
      const RealT inv_tap = ONE<RealT> / tap_;
      const RealT cos_ph  = std::cos(phase_);
      const RealT sin_ph  = std::sin(phase_);

      g11_ = -g_br * inv_tap * inv_tap - HALF<RealT> * G_ - Gmag_;
      b11_ = -b_br * inv_tap * inv_tap - HALF<RealT> * B_ - Bmag_;

      g12_ = (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = (b_br * cos_ph + g_br * sin_ph) * inv_tap;

      g21_ = (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = (b_br * cos_ph - g_br * sin_ph) * inv_tap;

      g22_ = -g_br - HALF<RealT> * G_;
      b22_ = -b_br - HALF<RealT> * B_;

      // An out-of-service branch has zero admittance
      g11_ *= in_service_;
      b11_ *= in_service_;
      g12_ *= in_service_;
      b12_ *= in_service_;
      g21_ *= in_service_;
      b21_ *= in_service_;
      g22_ *= in_service_;
      b22_ *= in_service_;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
