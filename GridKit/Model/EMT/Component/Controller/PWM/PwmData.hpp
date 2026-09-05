#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      enum class PwmParameters
      {
        M,
        fm,
        fc,
        alignment,
      };

      enum class PwmInputs : size_t
      {
        SIZE,
      };

      enum class PwmOutputs : size_t
      {
        sa,
        sb,
        sc,
        SIZE,
      };

      enum class PwmMonitorableVariables
      {
        sa,
        sb,
        sc,
      };

      template <typename real_type, typename index_type>
      struct PwmData : public ComponentData<real_type,
                                            index_type,
                                            PwmParameters,
                                            PwmInputs,
                                            PwmOutputs,
                                            PwmMonitorableVariables>
      {
        PwmData() = default;

        using Parameters           = PwmParameters;
        using Inputs               = PwmInputs;
        using Outputs              = PwmOutputs;
        using MonitorableVariables = PwmMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
