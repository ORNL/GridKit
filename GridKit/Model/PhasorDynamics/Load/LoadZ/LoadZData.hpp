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

    /// Buses for a load
    enum class LoadZBuses : size_t
    {
      bus, ///< Unique ID of the bus to which the load is connected
      SIZE
    };

    /// Signal inputs supported for a load
    enum class LoadZSignalInputs : size_t
    {
      SIZE
    };

    /// Signal outputs supported for a load
    enum class LoadZSignalOutputs : size_t
    {
      SIZE
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
    using LoadZData =
        ComponentData<real_type,
                      index_type,
                      LoadZParameters,
                      LoadZBuses,
                      LoadZSignalInputs,
                      LoadZSignalOutputs,
                      LoadZMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
