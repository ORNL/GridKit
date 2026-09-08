#pragma once

#include <GridKit/Model/EMT/Operators/Converter/ConverterData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const json& raw, ConverterData<RealT, IdxT>& data)
    {
      auto j = raw;
      expandPhasePort(j, "outputs", "e", {"ea", "eb", "ec"});
      expandPhasePort(j, "inputs", "s", {"sa", "sb", "sc"});
      expandPhasePort(j, "inputs", "i", {"ia", "ib", "ic"});
      expandPhaseMonitor(j, "e", {"ea", "eb", "ec"});
      using BaseT = ComponentData<RealT, IdxT, ConverterParameters, ConverterInputs, ConverterOutputs, ConverterMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
