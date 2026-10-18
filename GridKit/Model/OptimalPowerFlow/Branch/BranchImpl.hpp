/**
 * @file BranchImpl.hpp
 * @brief Definition of an optimal power flow branch.
 */

#pragma once

#include <array>
#include <cmath>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::Branch(const ModelDataT& data)
      : ComponentT(std::type_identity<Branch>{}, data.id)
    {
      initializeParameters(data);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Branch<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Branch");
      checks.check(R_ * R_ + X_ * X_ > ZERO<RealT>, "R and X cannot both be zero");
      checks.check(tap_ > ZERO<RealT>, "tap must be positive");
      checks.check(Smax_ > ZERO<RealT>, "Smax must be positive");
      return checks;
    }

    /**
     * @brief The squared apparent power at each port is an inequality
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::tagInequality()
    {
      const auto S1 = static_cast<size_t>(BranchInternalConstraints::S1);
      const auto S2 = static_cast<size_t>(BranchInternalConstraints::S2);

      tag_[S1] = true;
      tag_[S2] = true;
      return 0;
    }

    /**
     * @brief The rating bounds the squared apparent power at each port
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::initialize()
    {
      const auto S1 = static_cast<size_t>(BranchInternalConstraints::S1);
      const auto S2 = static_cast<size_t>(BranchInternalConstraints::S2);

      RealT* g_upper = g_upper_.getData();

      g_upper[S1] = Smax_ * Smax_;
      g_upper[S2] = Smax_ * Smax_;
      return 0;
    }

    /**
     * @brief Squared apparent power at bus 1 and at bus 2
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Branch<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                         const ScalarT*                  external,
                                                         ScalarT*                        g) const
    {
      const auto DP1 = static_cast<size_t>(BranchExternalConstraints::DP1);
      const auto DQ1 = static_cast<size_t>(BranchExternalConstraints::DQ1);
      const auto DP2 = static_cast<size_t>(BranchExternalConstraints::DP2);
      const auto DQ2 = static_cast<size_t>(BranchExternalConstraints::DQ2);
      const auto S1  = static_cast<size_t>(BranchInternalConstraints::S1);
      const auto S2  = static_cast<size_t>(BranchInternalConstraints::S2);

      std::array<ScalarT, Utilities::enum_size<BranchExternalConstraints>()> s;
      power(external, s.data());

      g[S1] = s[DP1] * s[DP1] + s[DQ1] * s[DQ1];
      g[S2] = s[DP2] * s[DP2] + s[DQ2] * s[DQ2];
    }

    /**
     * @brief Active and reactive power into bus 1, then into bus 2, added to
     * their balances
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Branch<scalar_type, index_type>::externalConstraints([[maybe_unused]] const ScalarT* x,
                                                         const ScalarT*                  external,
                                                         ScalarT*                        g) const
    {
      power(external, g);
    }

    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = BranchParameters;

      const Model::ParameterReader reader(data, "Branch");

      R_          = reader.readReal(Params::R).value_or(R_);
      X_          = reader.readReal(Params::X).value_or(X_);
      G_          = reader.readReal(Params::G).value_or(G_);
      B_          = reader.readReal(Params::B).value_or(B_);
      Gmag_       = reader.readReal(Params::Gmag).value_or(Gmag_);
      Bmag_       = reader.readReal(Params::Bmag).value_or(Bmag_);
      tap_        = reader.readReal(Params::tap).value_or(tap_);
      phase_      = reader.readReal(Params::phase).value_or(phase_);
      Smax_       = reader.readReal(Params::Smax).value_or(Smax_);
      in_service_ = reader.readBool(Params::in_service).value_or(in_service_);
      setDerivedParameters();
    }

    /**
     * @brief Pi-model admittance with the off-nominal tap and phase shift on
     * the bus-1 side, zero out of service
     */
    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::setDerivedParameters()
    {
      RealT u = ZERO<RealT>;
      if (in_service_)
      {
        u = ONE<RealT>;
      }

      const RealT denom   = R_ * R_ + X_ * X_;
      const RealT g_br    = R_ / denom;
      const RealT b_br    = -X_ / denom;
      const RealT inv_tap = ONE<RealT> / tap_;
      const RealT cos_ph  = std::cos(phase_);
      const RealT sin_ph  = std::sin(phase_);

      g11_ = u * (-g_br * inv_tap * inv_tap - HALF<RealT> * G_ - Gmag_);
      b11_ = u * (-b_br * inv_tap * inv_tap - HALF<RealT> * B_ - Bmag_);
      g12_ = u * (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = u * (b_br * cos_ph + g_br * sin_ph) * inv_tap;
      g21_ = u * (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = u * (b_br * cos_ph - g_br * sin_ph) * inv_tap;
      g22_ = u * (-g_br - HALF<RealT> * G_);
      b22_ = u * (-b_br - HALF<RealT> * B_);
    }

    /**
     * @brief Active and reactive power into bus 1, then into bus 2
     *
     * The power into each bus is \f$P + jQ = V I^*\f$.
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Branch<scalar_type, index_type>::power(const ScalarT* external, ScalarT* s) const
    {
      const auto VR1 = static_cast<size_t>(BranchExternalVariables::VR1);
      const auto VI1 = static_cast<size_t>(BranchExternalVariables::VI1);
      const auto VR2 = static_cast<size_t>(BranchExternalVariables::VR2);
      const auto VI2 = static_cast<size_t>(BranchExternalVariables::VI2);
      const auto DP1 = static_cast<size_t>(BranchExternalConstraints::DP1);
      const auto DQ1 = static_cast<size_t>(BranchExternalConstraints::DQ1);
      const auto DP2 = static_cast<size_t>(BranchExternalConstraints::DP2);
      const auto DQ2 = static_cast<size_t>(BranchExternalConstraints::DQ2);

      const ScalarT vr1 = external[VR1];
      const ScalarT vi1 = external[VI1];
      const ScalarT vr2 = external[VR2];
      const ScalarT vi2 = external[VI2];

      const ScalarT ir1 = g11_ * vr1 - b11_ * vi1 + g12_ * vr2 - b12_ * vi2;
      const ScalarT ii1 = b11_ * vr1 + g11_ * vi1 + b12_ * vr2 + g12_ * vi2;
      const ScalarT ir2 = g21_ * vr1 - b21_ * vi1 + g22_ * vr2 - b22_ * vi2;
      const ScalarT ii2 = b21_ * vr1 + g21_ * vi1 + b22_ * vr2 + g22_ * vi2;

      s[DP1] = vr1 * ir1 + vi1 * ii1;
      s[DQ1] = vi1 * ir1 - vr1 * ii1;
      s[DP2] = vr2 * ir2 + vi2 * ii2;
      s[DQ2] = vi2 * ir2 - vr2 * ii2;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
