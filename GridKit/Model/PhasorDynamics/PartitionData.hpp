#pragma once

#include <cstddef>
#include <filesystem>
#include <istream>
#include <string>
#include <vector>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// A named set of buses that one system model solves for
    template <typename real_type = double, typename index_type = size_t>
    struct PartitionData
    {
      using IdxT = index_type;

      std::string       name;  ///< A name given to this partition
      std::vector<IdxT> buses; ///< IDs of the buses in this partition
    };

    std::vector<PartitionData<double, size_t>> parsePartitionData(std::istream& stream);
    std::vector<PartitionData<double, size_t>> parsePartitionData(const std::filesystem::path& file);
  } // namespace PhasorDynamics
} // namespace GridKit
