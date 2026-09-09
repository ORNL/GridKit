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
        expandPhasePort<2>(j, "inputs", "u", {"ud", "uq"});
        expandPhasePort(j, "outputs", "s", {"sa", "sb", "sc"});
        expandPhasePort<2>(j, "outputs", "ulim", {"ulimd", "ulimq"});
        expandPhaseMonitor(j, "s", {"sa", "sb", "sc"});
        expandPhaseMonitor(j, "m", {"ma", "mb", "mc"});
        expandPhaseMonitor<2>(j, "ulim", {"ulimd", "ulimq"});
        using BaseT = ComponentData<RealT, IdxT, PwmParameters, PwmInputs, PwmOutputs, PwmMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
