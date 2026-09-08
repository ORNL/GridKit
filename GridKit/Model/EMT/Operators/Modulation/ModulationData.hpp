#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ModulationParameters
    {
    };

    enum class ModulationInputs : size_t
    {
      ua,
      ub,
      uc,
      vdc,
      SIZE,
    };

    enum class ModulationOutputs : size_t
    {
      ma,
      mb,
      mc,
      SIZE,
    };

    enum class ModulationMonitorableVariables
    {
      ma,
      mb,
      mc,
    };

    template <typename real_type, typename index_type>
    struct ModulationData : public ComponentData<real_type,
                                                 index_type,
                                                 ModulationParameters,
                                                 ModulationInputs,
                                                 ModulationOutputs,
                                                 ModulationMonitorableVariables>
    {
      ModulationData() = default;

      using Parameters           = ModulationParameters;
      using Inputs               = ModulationInputs;
      using Outputs              = ModulationOutputs;
      using MonitorableVariables = ModulationMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
