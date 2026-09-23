/**
 * @file ShuntData.hpp
 * @brief Modeling data for optimal power flow shunts.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a shunt
    enum class ShuntParameters : size_t
    {
      G, ///< \f$G\f$ Shunt conductance [p.u.]
      B, ///< \f$B\f$ Shunt susceptance [p.u.]
    };

    /// Buses for a shunt
    enum class ShuntBuses : size_t
    {
      bus, ///< Terminal bus
    };

    /**
     * @brief Contains modeling data for a Shunt
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    using ShuntData = ComponentData<real_type, index_type, ShuntParameters, ShuntBuses>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
