#pragma once

#include <GridKit/Model/EMT/Operators/Reference/Park/ParkData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const json& raw, ParkData<RealT, IdxT>& data)
    {
      auto j = raw;
      expandPhasePort(j, "inputs", "input", {"u1", "u2", "u3"});
      expandPhasePort(j, "outputs", "out", {"y1", "y2", "y3"});
      expandPhaseMonitor(j, "out", {"y1", "y2", "y3"});
      using BaseT = ComponentData<RealT, IdxT, ParkParameters, ParkInputs, ParkOutputs, ParkMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
