#pragma once

#include <GridKit/Model/EMT/Component/Line/LineDistributed/LineDistributedData.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Propagation/PropagationDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const nlohmann::json& j, LineDistributedData<RealT, IdxT>& data)
    {
      using BaseT = ComponentData<RealT, IdxT, LineDistributedParameters, LineDistributedInputs, LineDistributedOutputs, LineDistributedMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data), {"submodels"});
      const auto& submodels = j.at("submodels");
      validateJsonFields(submodels, "LineDistributed submodels", {"Yc", "H"});
      data.Yc = parseVectorFitOperand<RealT, IdxT>(submodels.at("Yc"));
      data.H  = parsePropagationOperand<RealT, IdxT>(submodels.at("H"));
    }
  } // namespace EMT
} // namespace GridKit
