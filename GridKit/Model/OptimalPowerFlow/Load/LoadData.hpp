/**
 * @file LoadData.hpp
 * @brief Modeling data for optimal power flow loads.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a load, whose demand comes from the state
    enum class LoadParameters : size_t
    {
    };

    /// Buses for a load
    enum class LoadBuses : size_t
    {
      bus, ///< Terminal bus
    };

    /**
     * @brief Contains modeling data for a Load
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    using LoadData = ComponentData<real_type, index_type, LoadParameters, LoadBuses>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
