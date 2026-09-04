#pragma once

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSourceData.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// JSON parser implementation for `DependentVoltageSourceData`
    ///
    /// The `submodels` object carries the rational source admittance under
    /// `Y`, as a sidecar fit file path or an inline object.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, DependentVoltageSourceData<RealT, IdxT>& d)
    {
      using BaseT = ComponentData<RealT,
                                  IdxT,
                                  DependentVoltageSourceParameters,
                                  DependentVoltageSourceBuses,
                                  DependentVoltageSourceSignalInputs,
                                  DependentVoltageSourceSignalOutputs,
                                  DependentVoltageSourceMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(d));

      if (j.contains("submodels") && j.at("submodels").contains("Y"))
      {
        d.Y = parseVectorFitOperand<RealT, IdxT>(j.at("submodels").at("Y"));
      }
    }
  } // namespace EMT
} // namespace GridKit
