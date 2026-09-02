#include "SystemModelData.hpp"

#include <fstream>

#include <GridKit/Model/EMT/SystemModelDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
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
      if (!stream)
      {
        std::stringstream ss;
        ss << "Could not open file: " << filePath;
        Log::error() << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      return parseSystemModelData(stream);
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::string& fileName)
    {
      auto stream = std::ifstream(fileName);
      if (!stream)
      {
        std::stringstream ss;
        ss << "Could not open file: " << fileName;
        Log::error() << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      return parseSystemModelData(stream);
    }
  } // namespace EMT
} // namespace GridKit
