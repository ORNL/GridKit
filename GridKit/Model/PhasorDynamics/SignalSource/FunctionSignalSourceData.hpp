#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    enum class FunctionSignalSourceParameters : size_t
    {
      Fr,
      Fi
    };

    enum class FunctionSignalSourceBuses : size_t
    {
    };

    enum class FunctionSignalSourceSignalInputs : size_t
    {
    };

    enum class FunctionSignalSourceSignalOutputs : size_t
    {
      sr,
      si
    };

    enum class FunctionSignalSourceMonitorableVariables : size_t
    {
    };

    template <typename real_type, typename index_type>
    using FunctionSignalSourceData =
        ComponentData<real_type,
                      index_type,
                      FunctionSignalSourceParameters,
                      FunctionSignalSourceBuses,
                      FunctionSignalSourceSignalInputs,
                      FunctionSignalSourceSignalOutputs,
                      FunctionSignalSourceMonitorableVariables>;

  } // namespace PhasorDynamics
} // namespace GridKit
