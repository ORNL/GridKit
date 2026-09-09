#pragma once

#include <GridKit/Model/EMT/Component/Controller/REPCA/RepcaData.hpp>
#include <GridKit/Model/EMT/PhaseSignalDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename RealT, typename IdxT>
      void from_json(const json& raw, RepcaData<RealT, IdxT>& data)
      {
        auto j = raw;
        expandPhasePort<2>(j, "inputs", "v", {"vd", "vq"});
        expandPhasePort<2>(j, "inputs", "i", {"id", "iq"});
        using BaseT = ComponentData<RealT, IdxT, RepcaParameters, RepcaInputs, RepcaOutputs, RepcaMonitorableVariables>;
        from_json(j, static_cast<BaseT&>(data));
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
