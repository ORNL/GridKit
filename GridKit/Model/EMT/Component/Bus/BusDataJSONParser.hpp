#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>

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

      if (j.contains("init"))
      {
        const auto& init = j.at("init");
        data.va0         = init.value("va", RealT{0});
        data.vb0         = init.value("vb", RealT{0});
        data.vc0         = init.value("vc", RealT{0});
      }
    }
  } // namespace EMT
} // namespace GridKit
