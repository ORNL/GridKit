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
      expandPhasePort(j, "inputs", "u", {"ua", "ub", "uc"});
      expandPhasePort(j, "outputs", "m", {"ma", "mb", "mc"});
      expandPhaseMonitor(j, "m", {"ma", "mb", "mc"});
      using BaseT = ComponentData<RealT, IdxT, ModulationParameters, ModulationInputs, ModulationOutputs, ModulationMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
