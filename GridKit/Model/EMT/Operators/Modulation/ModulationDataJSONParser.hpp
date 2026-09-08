#pragma once

#include <GridKit/Model/EMT/Operators/Modulation/ModulationData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const json& raw, ModulationData<RealT, IdxT>& data)
    {
      auto j = raw;
      expandPhasePort<2>(j, "inputs", "u", {"ud", "uq"});
      expandPhasePort<2>(j, "outputs", "m", {"md", "mq"});
      expandPhasePort<2>(j, "outputs", "ulim", {"ulimd", "ulimq"});
      expandPhaseMonitor<2>(j, "m", {"md", "mq"});
      expandPhaseMonitor<2>(j, "ulim", {"ulimd", "ulimq"});
      using BaseT = ComponentData<RealT, IdxT, ModulationParameters, ModulationInputs, ModulationOutputs, ModulationMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
