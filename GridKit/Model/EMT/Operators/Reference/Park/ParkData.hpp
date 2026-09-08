#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ParkParameters
    {
      inverse,
    };

    enum class ParkInputs : size_t
    {
      u1,
      u2,
      u3,
      theta,
      SIZE,
    };

    enum class ParkOutputs : size_t
    {
      y1,
      y2,
      y3,
      SIZE,
    };

    enum class ParkMonitorableVariables
    {
      y1,
      y2,
      y3,
    };

    template <typename real_type, typename index_type>
    struct ParkData : public ComponentData<real_type,
                                           index_type,
                                           ParkParameters,
                                           ParkInputs,
                                           ParkOutputs,
                                           ParkMonitorableVariables>
    {
      ParkData() = default;

      using Parameters           = ParkParameters;
      using Inputs               = ParkInputs;
      using Outputs              = ParkOutputs;
      using MonitorableVariables = ParkMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
