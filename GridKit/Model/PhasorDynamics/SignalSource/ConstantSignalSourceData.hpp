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

    enum class ConstantSignalSourceBuses : size_t
    {
      SIZE
    };

    enum class ConstantSignalSourceSignalInputs : size_t
    {
      SIZE
    };

    enum class ConstantSignalSourceSignalOutputs : size_t
    {
      sr,
      si,
      SIZE
    };

    enum class ConstantSignalSourceMonitorableVariables
    {
      NONE,
    };

    template <typename real_type, typename index_type>
    struct ConstantSignalSourceData : public ComponentData<real_type,
                                                           index_type,
                                                           ConstantSignalSourceParameters,
                                                           ConstantSignalSourceBuses,
                                                           ConstantSignalSourceSignalInputs,
                                                           ConstantSignalSourceSignalOutputs,
                                                           ConstantSignalSourceMonitorableVariables>
    {
      ConstantSignalSourceData() = default;

      using Parameters           = ConstantSignalSourceParameters;
      using Buses                = ConstantSignalSourceBuses;
      using SignalInputs         = ConstantSignalSourceSignalInputs;
      using SignalOutputs        = ConstantSignalSourceSignalOutputs;
      using MonitorableVariables = ConstantSignalSourceMonitorableVariables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
