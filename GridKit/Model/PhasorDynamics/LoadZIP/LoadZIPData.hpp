
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a loadZIP
    enum class LoadZIPParameters
    {
      P0, ///< Load nominal real power
      Q0, ///< Load nominal reactive power
      V0, ///< Load nominal reactive power
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
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename RealT, typename IdxT>
    struct LoadZIPData : public ComponentData<RealT,
                                           IdxT,
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
