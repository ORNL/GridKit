#pragma once

#include <GridKit/Model/EMT/Component/Controller/PWM/PwmData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename RealT, typename IdxT>
      void from_json(const json& raw, PwmData<RealT, IdxT>& data)
      {
        auto j = raw;
        expandPhasePort(j, "inputs", "m", {"ma", "mb", "mc"});
        expandPhasePort(j, "outputs", "s", {"sa", "sb", "sc"});
        expandPhaseMonitor(j, "s", {"sa", "sb", "sc"});
        using BaseT = ComponentData<RealT, IdxT, PwmParameters, PwmInputs, PwmOutputs, PwmMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
