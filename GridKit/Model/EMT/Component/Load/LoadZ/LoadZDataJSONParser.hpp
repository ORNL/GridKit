#pragma once

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// JSON parser function implementation for the `LoadZData` type
    ///
    /// The `submodels` object carries the rational impedance under `Z`, as a
    /// sidecar fit file path or an inline object.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, LoadZData<RealT, IdxT>& d)
    {
      using BaseT = ComponentData<RealT,
                                  IdxT,
                                  LoadZParameters,
                                  LoadZBuses,
                                  LoadZSignalInputs,
                                  LoadZSignalOutputs,
                                  LoadZMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(d));

      if (j.contains("submodels") && j.at("submodels").contains("Z"))
      {
        d.Z = parseVectorFitOperand<RealT, IdxT>(j.at("submodels").at("Z"));
      }
    }
  } // namespace EMT
} // namespace GridKit
