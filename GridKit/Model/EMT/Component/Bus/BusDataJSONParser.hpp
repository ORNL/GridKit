#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// JSON parser function implementation for the `BusData` type
    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& data)
    {
      using BaseT = ComponentData<RealT,
                                  IdxT,
                                  BusParameters,
                                  BusInputs,
                                  BusOutputs,
                                  BusMonitorableVariables>;
      from_json(j, static_cast<BaseT&>(data));

      if (data.device_class != "Bus")
      {
        throw std::runtime_error("JSON parser failed: expected Bus class");
      }
      data.shunts.clear();
      if (j.contains("shunts"))
      {
        if (!j.at("shunts").is_object())
          throw std::invalid_argument("Bus shunts must be a named object");
        for (const auto& [name, coefficients] : j.at("shunts").items())
          data.shunts.emplace(name, parseVectorFitOperand<RealT, IdxT>(coefficients));
      }
    }
  } // namespace EMT
} // namespace GridKit
