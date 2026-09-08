#pragma once

#include <GridKit/Model/EMT/Component/Source/REGFMA/RegfmaData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const json& raw, RegfmaData<RealT, IdxT>& data)
    {
      auto j = raw;
      expandPhasePort(j, "inputs", "v", {"va", "vb", "vc"});
      expandPhasePort(j, "outputs", "i", {"ia", "ib", "ic"});
      expandPhaseMonitor(j, "i", {"ia", "ib", "ic"});
      expandPhaseMonitor(j, "e", {"ea", "eb", "ec"});
      using BaseT = ComponentData<RealT, IdxT, RegfmaParameters, RegfmaInputs, RegfmaOutputs, RegfmaMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
