/**
 * @file LoadImpl.hpp
 * @brief Definition of an optimal power flow load.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentModelImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Load<scalar_type, index_type>::Load(const ModelDataT& data)
      : ComponentModelT(data.id, busNumbers(data))
    {
    }

    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::verify() const
    {
      return 0;
    }

    /**
     * @brief Take the demand from the state injection
     *
     * An offline load draws no power.
     */
    template <typename scalar_type, typename index_type>
    int Load<scalar_type, index_type>::initialize(const Model::StateData&   state,
                                                  [[maybe_unused]] ScalarT* x,
                                                  [[maybe_unused]] RealT*   x_lower,
                                                  [[maybe_unused]] RealT*   x_upper)
    {
      const bool has_power = this->statePower(state, 0, p_, q_);

      const Model::StateRecord* record = state.device(id_);
      if (record != nullptr && !record->flag("online", true))
      {
        p_ = ZERO<RealT>;
        q_ = ZERO<RealT>;
        return 0;
      }

      return this->check(has_power, "state has no current for the load or no voltage for its bus");
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
      g[0] = p_;
      g[1] = q_;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
