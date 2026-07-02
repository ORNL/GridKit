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

    enum class ConstantSignalSourceTerminals : size_t
    {
      SIZE
    };

    enum class ConstantSignalSourceInputPorts : size_t
    {
      SIZE
    };

    enum class ConstantSignalSourceOutputPorts : size_t
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
                                                           ConstantSignalSourceTerminals,
                                                           ConstantSignalSourceInputPorts,
                                                           ConstantSignalSourceOutputPorts,
                                                           ConstantSignalSourceMonitorableVariables>
    {
      ConstantSignalSourceData() = default;

      using Parameters           = ConstantSignalSourceParameters;
      using Terminals            = ConstantSignalSourceTerminals;
      using InputPorts           = ConstantSignalSourceInputPorts;
      using OutputPorts          = ConstantSignalSourceOutputPorts;
      using MonitorableVariables = ConstantSignalSourceMonitorableVariables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
