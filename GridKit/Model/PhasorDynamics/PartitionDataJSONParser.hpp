#pragma once

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/PartitionData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `PartitionData` type
    template <typename RealT, typename IdxT>
    void from_json(const json& j, PartitionData<RealT, IdxT>& pd)
    {
      j.at("name").get_to(pd.name);
      j.at("buses").get_to(pd.buses);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
