/**
 * @file LoadImpl.hpp
 * @brief Definition of an optimal power flow load.
 */

#pragma once

#include <cmath>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(const ModelDataT& data)
      : ComponentT(std::type_identity<Load>{}, data.id)
    {
      initializeParameters(data);
    }

    /// Every demand is valid
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Load<scalar_type, index_type>::verify() const
    {
      return Model::ConfigurationChecks("Load");
    }

    /// Power into the bus, the negative of the demand
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Load<scalar_type, index_type>::externalConstraints([[maybe_unused]] const ScalarT* x,
                                                       const ScalarT*                  external,
                                                       ScalarT*                        g) const
    {
      const auto VR = static_cast<size_t>(LoadExternalVariables::VR);
      const auto VI = static_cast<size_t>(LoadExternalVariables::VI);
      const auto DP = static_cast<size_t>(LoadExternalConstraints::DP);
      const auto DQ = static_cast<size_t>(LoadExternalConstraints::DQ);

      const ScalarT v2 = external[VR] * external[VR] + external[VI] * external[VI];
      const ScalarT v  = std::sqrt(current_on_ * v2 + (ONE<RealT> - current_on_));

      g[DP] = -u_ * (P_ + Pi_ * v + Pz_ * v2);
      g[DQ] = -u_ * (Q_ + Qi_ * v + Qz_ * v2);
    }

    template <typename scalar_type, typename index_type>
    void Load<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = LoadParameters;

      const Model::ParameterReader reader(data, "Load");

      P_          = reader.readReal(Params::P).value_or(P_);
      Q_          = reader.readReal(Params::Q).value_or(Q_);
      Pi_         = reader.readReal(Params::Pi).value_or(Pi_);
      Qi_         = reader.readReal(Params::Qi).value_or(Qi_);
      Pz_         = reader.readReal(Params::Pz).value_or(Pz_);
      Qz_         = reader.readReal(Params::Qz).value_or(Qz_);
      in_service_ = reader.readBool(Params::in_service).value_or(in_service_);
      setDerivedParameters();
    }

    template <typename scalar_type, typename index_type>
    void Load<scalar_type, index_type>::setDerivedParameters()
    {
      u_ = ZERO<RealT>;
      if (in_service_)
      {
        u_ = ONE<RealT>;
      }

      const bool current_part = Pi_ < ZERO<RealT> || Pi_ > ZERO<RealT> || Qi_ < ZERO<RealT> || Qi_ > ZERO<RealT>;

      current_on_ = ZERO<RealT>;
      if (in_service_ && current_part)
      {
        current_on_ = ONE<RealT>;
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
