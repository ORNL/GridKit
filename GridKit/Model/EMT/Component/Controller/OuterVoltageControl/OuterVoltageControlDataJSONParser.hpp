#pragma once

#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControlData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename RealT, typename IdxT>
      void from_json(const json& raw, OuterVoltageControlData<RealT, IdxT>& data)
      {
        auto j = raw;
        expandPhasePort<2>(j, "inputs", "vref", {"vrefd", "vrefq"});
        expandPhasePort<2>(j, "inputs", "v", {"vd", "vq"});
        expandPhasePort<2>(j, "inputs", "ig", {"igd", "igq"});
        expandPhasePort<2>(j, "inputs", "ilim", {"ilimd", "ilimq"});
        expandPhasePort<2>(j, "outputs", "iref", {"irefd", "irefq"});
        expandPhaseMonitor<2>(j, "eta", {"etad", "etaq"});
        expandPhaseMonitor<2>(j, "iref", {"irefd", "irefq"});
        using BaseT = ComponentData<RealT, IdxT, OuterVoltageControlParameters, OuterVoltageControlInputs, OuterVoltageControlOutputs, OuterVoltageControlMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
