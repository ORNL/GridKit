#include "SystemModelData.hpp"

#include <fstream>

#include <Model/PhasorDynamics/SystemModelDataJSONParser.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream)
    {
      SystemModelData<double, size_t> data(json::parse(stream));
      return data;
    }

    SystemModelData<double, size_t> parseSystemModelData(std::istream&& stream)
    {
      SystemModelData<double, size_t> data(json::parse(stream));
      return data;
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& filePath)
    {
      auto stream = std::ifstream(filePath);
      return parseSystemModelData(stream);
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::string& fileName)
    {
      return parseSystemModelData(std::ifstream(fileName));
    }
  } // namespace PhasorDynamics
} // namespace GridKit
