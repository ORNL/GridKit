#include "PartitionData.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include <GridKit/Model/PhasorDynamics/PartitionDataJSONParser.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    std::vector<PartitionData<double, size_t>> parsePartitionData(std::istream& stream)
    {
      return json::parse(stream).at("partitions").get<std::vector<PartitionData<double, size_t>>>();
    }

    std::vector<PartitionData<double, size_t>> parsePartitionData(const std::filesystem::path& file)
    {
      auto stream = std::ifstream(file);
      if (!stream)
      {
        std::stringstream ss;
        ss << "Could not open file: " << file;
        Log::error() << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      return parsePartitionData(stream);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
