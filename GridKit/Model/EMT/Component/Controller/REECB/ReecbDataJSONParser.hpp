#pragma once

#include <GridKit/Model/EMT/Component/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename RealT, typename IdxT>
      void from_json(const json& raw, ReecbData<RealT, IdxT>& data)
      {
        auto j = raw;
        expandPhasePort<2>(j, "inputs", "v", {"vd", "vq"});
        expandPhasePort<2>(j, "inputs", "i", {"id", "iq"});
        expandPhasePort<2>(j, "outputs", "icmd", {"icmdd", "icmdq"});
        expandPhaseMonitor<2>(j, "icmd", {"icmdd", "icmdq"});
        using BaseT = ComponentData<RealT, IdxT, ReecbParameters, ReecbInputs, ReecbOutputs, ReecbMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
