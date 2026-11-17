#include "PartitionData.hpp"

#include <fstream>

#include <GridKit/Model/PhasorDynamics/Partition/PartitionDataJSONParser.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    PartitionData<size_t> parsePartitionData(std::istream& stream)
    {
      PartitionData<size_t> data;
      json::parse(stream).get_to(data);
      return data;
    }

    PartitionData<size_t> parsePartitionData(const std::filesystem::path& file_path)
    {
      auto stream = std::ifstream(file_path);
      if (!stream)
      {
        Log::error() << "Could not open partition file: " << file_path << "\n";
        throw std::runtime_error("JSON parser failed");
      }
      return parsePartitionData(stream);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
