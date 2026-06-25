/**
 * @file LoadZData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for loads
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a load
    enum class LoadZParameters
    {
      R, ///< Load resistance
      X, ///< Load reactance
    };

    /// Ports for a load
    enum class LoadZPorts
    {
      bus, ///< Unique ID of the bus to which the load is connected
    };

    /// Variables able to be monitored for a load
    enum class LoadZMonitorableVariables
    {
      p,
      q
    };

    /**
     * @brief Contains modeling data for a load
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct LoadZData : public ComponentData<real_type,
                                            index_type,
                                            LoadZParameters,
                                            LoadZPorts,
                                            LoadZMonitorableVariables>
    {
      LoadZData() = default;

      using Parameters           = LoadZParameters;
      using Ports                = LoadZPorts;
      using MonitorableVariables = LoadZMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
