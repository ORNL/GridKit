#pragma once

#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControlData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename RealT, typename IdxT>
      void from_json(const json& raw, InnerCurrentControlData<RealT, IdxT>& data)
      {
        auto j = raw;
        expandPhasePort<2>(j, "inputs", "v", {"vd", "vq"});
        expandPhasePort<2>(j, "inputs", "i", {"id", "iq"});
        expandPhasePort<2>(j, "inputs", "icmd", {"icmdd", "icmdq"});
        expandPhasePort<2>(j, "inputs", "ulim", {"ulimd", "ulimq"});
        expandPhasePort<2>(j, "outputs", "ilim", {"ilimd", "ilimq"});
        expandPhasePort<2>(j, "outputs", "u", {"ud", "uq"});
        expandPhaseMonitor<2>(j, "xi", {"xid", "xiq"});
        expandPhaseMonitor<2>(j, "ilim", {"ilimd", "ilimq"});
        expandPhaseMonitor<2>(j, "u", {"ud", "uq"});
        using BaseT = ComponentData<RealT, IdxT, InnerCurrentControlParameters, InnerCurrentControlInputs, InnerCurrentControlOutputs, InnerCurrentControlMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
