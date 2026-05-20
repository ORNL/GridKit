
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a loadZIP
    enum class LoadZIPParameters
    {
      P0,     ///< Load nominal real power
      Q0,     ///< Load nominal reactive power
      V0,     ///< Load nominal reactive power
      alphaI, ///< Fraction of load to be represented as constant current
      alphaP, ///< Fraction of load to be represented as constant power
    };

    /// Ports for a loadZIP
    enum class LoadZIPPorts
    {
      bus, ///< Unique ID of the bus to which the loadZIP is connected
    };

    /// Variables able to be monitored for a loadZIP
    enum class LoadZIPMonitorableVariables
    {
      Ir,
      Ii,
      Im,
      P,
      Q
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
    struct LoadZIPData : public ComponentData<real_type,
                                              index_type,
                                              LoadZIPParameters,
                                              LoadZIPPorts,
                                              LoadZIPMonitorableVariables>
    {
      LoadZIPData() = default;

      using Parameters           = LoadZIPParameters;
      using Ports                = LoadZIPPorts;
      using MonitorableVariables = LoadZIPMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
