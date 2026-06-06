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

    /// Terminals for a load
    enum class LoadZTerminals : size_t
    {
      bus, ///< Unique ID of the bus to which the load is connected
      SIZE
    };

    /// Input ports supported for a load
    enum class LoadZInputPorts : size_t
    {
      SIZE
    };

    /// Output ports supported for a load
    enum class LoadZOutputPorts : size_t
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
    struct LoadZData : public ComponentData<real_type,
                                            index_type,
                                            LoadZParameters,
                                            LoadZTerminals,
                                            LoadZInputPorts,
                                            LoadZOutputPorts,
                                            LoadZMonitorableVariables>
    {
      LoadZData() = default;

      using Parameters           = LoadZParameters;
      using Terminals            = LoadZTerminals;
      using InputPorts           = LoadZInputPorts;
      using OutputPorts          = LoadZOutputPorts;
      using MonitorableVariables = LoadZMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
