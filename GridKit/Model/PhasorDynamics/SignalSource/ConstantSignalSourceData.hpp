#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    enum class ConstantSignalSourceParameters : size_t
    {
      Sr,
      Si,
    };

    enum class ConstantSignalSourceBuses : size_t
    {
    };

    enum class ConstantSignalSourceSignalInputs : size_t
    {
    };

    enum class ConstantSignalSourceSignalOutputs : size_t
    {
      sr,
      si,
    };

    enum class ConstantSignalSourceMonitorableVariables : size_t
    {
    };

    template <typename real_type, typename index_type>
    using ConstantSignalSourceData =
        ComponentData<real_type,
                      index_type,
                      ConstantSignalSourceParameters,
                      ConstantSignalSourceBuses,
                      ConstantSignalSourceSignalInputs,
                      ConstantSignalSourceSignalOutputs,
                      ConstantSignalSourceMonitorableVariables>;

  } // namespace PhasorDynamics
} // namespace GridKit
