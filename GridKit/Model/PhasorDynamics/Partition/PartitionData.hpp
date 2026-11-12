/**
 * @file PartitionData.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Contents of a `.partition.json` file.
 */
#pragma once

#include <filesystem>
#include <istream>
#include <string>
#include <vector>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief One partition of a partition file and its child partitions
     *
     * The file format is described in README.md in this directory.
     */
    template <typename index_type = size_t>
    struct PartitionData
    {
      using IdxT = index_type;

      std::string                name;
      std::vector<IdxT>          internal;   ///< Case bus IDs internal to this partition
      std::vector<IdxT>          external;   ///< Case bus IDs external to this partition
      std::vector<PartitionData> partitions; ///< Child partitions (none for a leaf)
    };

    /**
     * @brief Generate partition data from a JSON input stream or file
     *
     * @throws std::runtime_error for an unreadable file, an empty name or an
     * unknown key
     * @throws nlohmann::json exceptions for invalid JSON, a missing name or
     * wrong value types
     */
    PartitionData<size_t> parsePartitionData(std::istream& stream);
    PartitionData<size_t> parsePartitionData(const std::filesystem::path& file_path);
  } // namespace PhasorDynamics
} // namespace GridKit
