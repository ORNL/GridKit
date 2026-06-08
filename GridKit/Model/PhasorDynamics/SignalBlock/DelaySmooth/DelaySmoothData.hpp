#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalBlock
    {
      enum class DelaySmoothParameters
      {
        delay,  ///< Delay tau to approximate, in seconds (required, > 0)
        dt_min, ///< Block resolution; N = floor(tau/dt_min) (required, > 0)
      };

      enum class DelaySmoothPorts
      {
        input,
        output,
      };

      enum class DelaySmoothMonitorableVariables
      {
        out,
      };

      template <typename real_type, typename index_type>
      struct DelaySmoothData : public ComponentData<real_type,
                                                    index_type,
                                                    DelaySmoothParameters,
                                                    DelaySmoothPorts,
                                                    DelaySmoothMonitorableVariables>
      {
        using RealT                = real_type;
        using IdxT                 = index_type;
        using Parameters           = DelaySmoothParameters;
        using Ports                = DelaySmoothPorts;
        using MonitorableVariables = DelaySmoothMonitorableVariables;

        DelaySmoothData() = default;
      };
    } // namespace SignalBlock
  } // namespace PhasorDynamics
} // namespace GridKit
