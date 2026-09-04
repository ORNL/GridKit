#pragma once

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// JSON parser function implementation for the `LineLumpedData` type
    ///
    /// The `submodels` object carries the rational per-unit-length series
    /// impedance under `Zp` and the rational per-unit-length shunt admittance
    /// under `Yp`, each as a sidecar fit file path or an inline object.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, LineLumpedData<RealT, IdxT>& d)
    {
      using BaseT = ComponentData<RealT,
                                  IdxT,
                                  LineLumpedParameters,
                                  LineLumpedInputs,
                                  LineLumpedOutputs,
                                  LineLumpedMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(d));

      if (j.contains("submodels") && j.at("submodels").contains("Zp"))
      {
        d.Zp = parseVectorFitOperand<RealT, IdxT>(j.at("submodels").at("Zp"));
      }

      if (j.contains("submodels") && j.at("submodels").contains("Yp"))
      {
        d.Yp = parseVectorFitOperand<RealT, IdxT>(j.at("submodels").at("Yp"));
      }
    }
  } // namespace EMT
} // namespace GridKit
