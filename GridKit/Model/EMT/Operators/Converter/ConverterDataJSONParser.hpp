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
      expandPhasePort(j, "outputs", "vo", {"voa", "vob", "voc"});
      expandPhasePort(j, "inputs", "s", {"sa", "sb", "sc"});
      expandPhasePort(j, "inputs", "i", {"ia", "ib", "ic"});
      expandPhaseMonitor(j, "vo", {"voa", "vob", "voc"});
      using BaseT = ComponentData<RealT, IdxT, ConverterParameters, ConverterInputs, ConverterOutputs, ConverterMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
