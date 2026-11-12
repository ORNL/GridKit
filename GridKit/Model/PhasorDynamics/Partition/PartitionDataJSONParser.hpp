#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Partition/PartitionData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser function implementation for the `PartitionData` type
    ///
    /// Rejects empty names and unknown keys. See the `README.md` in
    /// `GridKit/Model/PhasorDynamics/Partition` for the file format.
    template <typename IdxT>
    void from_json(const json& j, PartitionData<IdxT>& pd)
    {
      j.at("name").get_to(pd.name);
      if (pd.name.empty())
      {
        Log::error() << "Partition names must not be empty\n";
        throw std::runtime_error("JSON parser failed");
      }

      for (const auto& entry : j.items())
      {
        if (entry.key() != "name" && entry.key() != "internal" && entry.key() != "external"
            && entry.key() != "partitions")
        {
          Log::error() << "Partition " << pd.name << " has unknown key " << entry.key() << "\n";
          throw std::runtime_error("JSON parser failed");
        }
      }

      if (j.contains("internal"))
      {
        j.at("internal").get_to(pd.internal);
      }

      if (j.contains("external"))
      {
        j.at("external").get_to(pd.external);
      }

      if (j.contains("partitions"))
      {
        j.at("partitions").get_to(pd.partitions);
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
