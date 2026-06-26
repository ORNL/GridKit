#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    enum class ConstantSignalSourceParameters
    {
      Sr,
      Si
    };

    enum class ConstantSignalSourcePorts
    {
      sr,
      si
    };

    enum class ConstantSignalSourceMonitorableVariables
    {
      NONE,
    };

    template <typename real_type, typename index_type>
    struct ConstantSignalSourceData : public ComponentData<real_type,
                                                           index_type,
                                                           ConstantSignalSourceParameters,
                                                           ConstantSignalSourcePorts,
                                                           ConstantSignalSourceMonitorableVariables>
    {
      ConstantSignalSourceData() = default;

      using Parameters           = ConstantSignalSourceParameters;
      using Ports                = ConstantSignalSourcePorts;
      using MonitorableVariables = ConstantSignalSourceMonitorableVariables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
