/**
 * @file LoadData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for loads
 *
 */
#pragma once

#include <Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a load
    enum class LoadParameters
    {
      /// Load resistance
      R,

      /// Load reactance
      X,
    };

    /// Ports for a load
    enum class LoadPorts
    {
      /// Unique ID of the bus to which the load is connected
      bus,
    };

    /// Variables able to be monitored for a load
    enum class LoadMonitorableVariables
    {
      // TODO: presumably some variables would make sense to monitor here
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
    struct LoadData : public ComponentData<RealT,
                                           IdxT,
                                           LoadParameters,
                                           LoadPorts,
                                           LoadMonitorableVariables>
    {
      LoadData() = default;

      using Parameters           = LoadParameters;
      using Ports                = LoadPorts;
      using MonitorableVariables = LoadMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
