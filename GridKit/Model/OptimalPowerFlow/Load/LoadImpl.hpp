/**
 * @file LoadImpl.hpp
 * @brief Definition of an optimal power flow load.
 */

#pragma once

#include <cmath>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(const ModelDataT& data)
      : ComponentT(std::type_identity<Load>{}, data.id),
        in_service_(data.in_service),
        current_part_(data.in_service
                      && (data.parameters.contains(LoadParameters::Pi) || data.parameters.contains(LoadParameters::Qi))),
        P_(parameter(data, LoadParameters::P, ZERO<RealT>)),
        Q_(parameter(data, LoadParameters::Q, ZERO<RealT>)),
        Pi_(parameter(data, LoadParameters::Pi, ZERO<RealT>)),
        Qi_(parameter(data, LoadParameters::Qi, ZERO<RealT>)),
        Pz_(parameter(data, LoadParameters::Pz, ZERO<RealT>)),
        Qz_(parameter(data, LoadParameters::Qz, ZERO<RealT>))
    {
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Load<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Load");
      checks.check(std::isfinite(P_), "P must be finite");
      checks.check(std::isfinite(Q_), "Q must be finite");
      checks.check(std::isfinite(Pi_), "Pi must be finite");
      checks.check(std::isfinite(Qi_), "Qi must be finite");
      checks.check(std::isfinite(Pz_), "Pz must be finite");
      checks.check(std::isfinite(Qz_), "Qz must be finite");
      return checks;
    }

    /**
     * @brief Active and reactive power into the bus, the negative of the demand
     *
     * A load without a constant-current part uses a unit magnitude in its
     * place, so it stays smooth at zero voltage.
     */
    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Load<scalar_type, index_type>::externalEqualities(const ScalarT* x, ScalarT* c) const
    {
      const ScalarT v2 = x[0] * x[0] + x[1] * x[1];
      const ScalarT v  = std::sqrt(current_part_ * v2 + (ONE<RealT> - current_part_));

      c[0] = -in_service_ * (P_ + Pi_ * v + Pz_ * v2);
      c[1] = -in_service_ * (Q_ + Qi_ * v + Qz_ * v2);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
