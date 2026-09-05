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
        if (j.contains("params"))
        {
          for (const char* name : {"M", "fm", "fc", "alignment"})
          {
            if (j.at("params").contains(name))
            {
              const auto& value = j.at("params").at(name);
              if (!value.is_number())
              {
                throw std::invalid_argument(std::string("PWM parameter ") + name + " must be numeric");
              }
              j["params"][name] = value.template get<RealT>();
            }
          }
        }
        expandPhasePort(j, "outputs", "s", {"sa", "sb", "sc"});
        expandPhaseMonitor(j, "s", {"sa", "sb", "sc"});
        using BaseT = ComponentData<RealT, IdxT, PwmParameters, PwmInputs, PwmOutputs, PwmMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
