#include "SystemModelData.hpp"

#include <fstream>
#include <sstream>

#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace
    {
      void resolveRelativePaths(SystemModelData<double, size_t>& data, const std::filesystem::path& base_path)
      {
        for (auto& source : data.sampled_signal_source)
        {
          if (!source.file.empty() && source.file.is_relative())
          {
            source.file = base_path / source.file;
          }
        }
      }
    } // namespace

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
      auto data = parseSystemModelData(stream);
      resolveRelativePaths(data, filePath.parent_path());
      return data;
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::string& fileName)
    {
      const auto filePath = std::filesystem::path(fileName);
      auto       stream   = std::ifstream(filePath);
      if (!stream)
      {
        std::stringstream ss;
        ss << "Could not open file: " << filePath;
        Log::error() << ss.str() << std::endl;
        throw std::runtime_error(ss.str());
      }
      auto data = parseSystemModelData(stream);
      resolveRelativePaths(data, filePath.parent_path());
      return data;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
