/**
 * @file LoadImpl.hpp
 * @brief Definition of an optimal power flow load.
 */

#pragma once

#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(const ModelDataT& data)
      : ComponentModelT(data.id, busNumbers(data), data.in_service),
        P_(parameter(data, LoadParameters::P, ZERO<RealT>)),
        Q_(parameter(data, LoadParameters::Q, ZERO<RealT>))
    {
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Load<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Load");
      checks.check(std::isfinite(P_), "P must be finite");
      checks.check(std::isfinite(Q_), "Q must be finite");
      return checks;
    }

    /**
     * @brief No internal variables
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::initialize([[maybe_unused]] ScalarT* x)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Load<scalar_type, index_type>::ScalarT
    Load<scalar_type, index_type>::objective([[maybe_unused]] const ScalarT* x) const
    {
      return ScalarT{ZERO<RealT>};
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Load<scalar_type, index_type>::constraints([[maybe_unused]] const ScalarT* x, ScalarT* g) const
    {
      g[0] = -in_service_ * P_;
      g[1] = -in_service_ * Q_;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
