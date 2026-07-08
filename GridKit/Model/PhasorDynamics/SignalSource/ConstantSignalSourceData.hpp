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
