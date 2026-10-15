/**
 * @file ShuntImpl.hpp
 * @brief Definition of an optimal power flow shunt.
 */

#pragma once

#include <cmath>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Shunt/Shunt.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Shunt<scalar_type, index_type>::Shunt(const ModelDataT& data)
      : ComponentModelT(data.id, busNumbers(data), data.in_service),
        G_(parameter(data, ShuntParameters::G, ZERO<RealT>)),
        B_(parameter(data, ShuntParameters::B, ZERO<RealT>))
    {
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Shunt<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Shunt");
      checks.check(std::isfinite(G_), "G must be finite");
      checks.check(std::isfinite(B_), "B must be finite");
      return checks;
    }

    /**
     * @brief No internal variables
     */
    template <typename scalar_type, typename index_type>
    int Shunt<scalar_type, index_type>::initialize([[maybe_unused]] ScalarT* x)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline typename Shunt<scalar_type, index_type>::ScalarT
    Shunt<scalar_type, index_type>::objective([[maybe_unused]] const ScalarT* x) const
    {
      return ScalarT{ZERO<RealT>};
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void
    Shunt<scalar_type, index_type>::constraints(const ScalarT* x, ScalarT* g) const
    {
      const ScalarT v2 = x[0] * x[0] + x[1] * x[1];

      g[0] = -in_service_ * G_ * v2;
      g[1] = in_service_ * B_ * v2;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
