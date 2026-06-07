#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      enum class DelayParameters
      {
        delay,      ///< Transport delay tau in seconds (required, > 0)
        prehistory, ///< Optional constant output for t < t0; default = input(t0)
      };

      enum class DelayPorts
      {
        input,
        output,
      };

      enum class DelayMonitorableVariables
      {
        out,
      };

      template <typename real_type, typename index_type>
      struct DelayData : public ComponentData<real_type,
                                              index_type,
                                              DelayParameters,
                                              DelayPorts,
                                              DelayMonitorableVariables>
      {
        using RealT                = real_type;
        using IdxT                 = index_type;
        using Parameters           = DelayParameters;
        using Ports                = DelayPorts;
        using MonitorableVariables = DelayMonitorableVariables;

        DelayData() = default;
      };
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
