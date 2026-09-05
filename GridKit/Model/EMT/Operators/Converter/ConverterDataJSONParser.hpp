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
      if (j.contains("params") && !j.at("params").empty())
      {
        throw std::invalid_argument("Converter has no parameters");
      }
      expandPhasePort(j, "outputs", "vo", {"voa", "vob", "voc"});
      expandPhasePort(j, "inputs", "s", {"sa", "sb", "sc"});
      expandPhaseMonitor(j, "vo", {"voa", "vob", "voc"});
      using BaseT = ComponentData<RealT, IdxT, ConverterParameters, ConverterInputs, ConverterOutputs, ConverterMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));
    }
  } // namespace EMT
} // namespace GridKit
