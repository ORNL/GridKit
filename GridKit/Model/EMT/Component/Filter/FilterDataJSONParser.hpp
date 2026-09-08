#pragma once

#include <GridKit/Model/EMT/Component/Filter/FilterData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const json& raw, FilterData<RealT, IdxT>& data)
    {
      auto j = raw;
      expandPhasePort(j, "inputs", "v", {"va", "vb", "vc"});
      expandPhasePort(j, "inputs", "e", {"ea", "eb", "ec"});
      expandPhasePort(j, "outputs", "i", {"ia", "ib", "ic"});
      expandPhasePort(j, "outputs", "vo", {"voa", "vob", "voc"});
      expandPhasePort(j, "outputs", "ig", {"iga", "igb", "igc"});
      expandPhaseMonitor(j, "i", {"ia", "ib", "ic"});
      expandPhaseMonitor(j, "vo", {"voa", "vob", "voc"});
      expandPhaseMonitor(j, "ig", {"iga", "igb", "igc"});
      using BaseT = ComponentData<RealT, IdxT, FilterParameters, FilterInputs, FilterOutputs, FilterMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
