/**
 * @file BranchImpl.hpp
 * @brief Definition of an optimal power flow branch.
 */

#pragma once

#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Construct a branch from data
     *
     * A terminal apparent power row exists only with a rating `Smax`.
     */
    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::Branch(const ModelDataT& data)
      : ComponentModelT(data.id, busNumbers(data)),
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
      constraint_upper_ = {Smax_ * Smax_, Smax_ * Smax_};
    }

    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::verify() const
    {
      return this->check(std::isfinite(R_), "R must be finite")
             + this->check(std::isfinite(X_), "X must be finite")
             + this->check(std::isfinite(G_), "G must be finite")
             + this->check(std::isfinite(B_), "B must be finite")
             + this->check(std::isfinite(Gmag_), "Gmag must be finite")
             + this->check(std::isfinite(Bmag_), "Bmag must be finite")
             + this->check(std::isfinite(tap_), "tap must be finite")
             + this->check(std::isfinite(phase_), "phase must be finite")
             + this->check(R_ * R_ + X_ * X_ > ZERO<RealT>, "R and X cannot both be zero")
             + this->check(tap_ > ZERO<RealT>, "tap must be positive")
             + this->check(Smax_ > ZERO<RealT>, "Smax must be positive");
    }

    /**
     * @brief Apply the state settings `open`, `tap` and `phase`
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::initialize(const Model::StateData&   state,
                                                    [[maybe_unused]] ScalarT* x,
                                                    [[maybe_unused]] RealT*   x_lower,
                                                    [[maybe_unused]] RealT*   x_upper)
    {
      RealT tap    = tap_;
      RealT phase  = phase_;
      RealT status = ONE<RealT>;

      const Model::StateRecord* record = state.device(id_);
      if (record != nullptr)
      {
        tap   = record->value("tap", tap);
        phase = record->value("phase", phase);
        if (record->flag("open", false))
        {
          status = ZERO<RealT>;
        }
      }

      setDerivedParameters(tap, phase, status);
      return this->check(tap > ZERO<RealT>, "state tap must be positive");
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Branch<scalar_type, index_type>::ScalarT
    Branch<scalar_type, index_type>::objective([[maybe_unused]] const ScalarT* x) const
    {
      return ScalarT{ZERO<RealT>};
    }

    /**
     * @brief Squared apparent power and power into each bus
     *
     * Terminal currents are those of `PhasorDynamics::Branch`, and the power
     * into bus k is \f$V_k I_k^*\f$.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Branch<scalar_type, index_type>::constraints(const ScalarT* x, ScalarT* g) const
    {
      const ScalarT vr1 = x[0];
      const ScalarT vi1 = x[1];
      const ScalarT vr2 = x[2];
      const ScalarT vi2 = x[3];

      const ScalarT ir1 = g11_ * vr1 - b11_ * vi1 + g12_ * vr2 - b12_ * vi2;
      const ScalarT ii1 = b11_ * vr1 + g11_ * vi1 + b12_ * vr2 + g12_ * vi2;
      const ScalarT ir2 = g21_ * vr1 - b21_ * vi1 + g22_ * vr2 - b22_ * vi2;
      const ScalarT ii2 = b21_ * vr1 + g21_ * vi1 + b22_ * vr2 + g22_ * vi2;

      const ScalarT p1 = vr1 * ir1 + vi1 * ii1;
      const ScalarT q1 = vi1 * ir1 - vr1 * ii1;
      const ScalarT p2 = vr2 * ir2 + vi2 * ii2;
      const ScalarT q2 = vi2 * ir2 - vr2 * ii2;

      g[0] = p1 * p1 + q1 * q1;
      g[1] = p2 * p2 + q2 * q2;
      g[2] = p1;
      g[3] = q1;
      g[4] = p2;
      g[5] = q2;
    }

    /**
     * @brief Admittance of `PhasorDynamics::Branch`, scaled by the branch status
     */
    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::setDerivedParameters(RealT tap, RealT phase, RealT status)
    {
      const RealT denom   = R_ * R_ + X_ * X_;
      const RealT g_br    = R_ / denom;
      const RealT b_br    = -X_ / denom;
      const RealT inv_tap = ONE<RealT> / tap;
      const RealT cos_ph  = std::cos(phase);
      const RealT sin_ph  = std::sin(phase);

      g11_ = status * (-g_br * inv_tap * inv_tap - HALF<RealT> * G_ - Gmag_);
      b11_ = status * (-b_br * inv_tap * inv_tap - HALF<RealT> * B_ - Bmag_);

      g12_ = status * (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = status * (b_br * cos_ph + g_br * sin_ph) * inv_tap;

      g21_ = status * (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = status * (b_br * cos_ph - g_br * sin_ph) * inv_tap;

      g22_ = status * (-g_br - HALF<RealT> * G_);
      b22_ = status * (-b_br - HALF<RealT> * B_);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
