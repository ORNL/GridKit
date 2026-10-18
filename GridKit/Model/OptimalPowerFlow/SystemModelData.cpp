/**
 * @file SystemModelData.cpp
 * @brief Optimal power flow data from a JSON input file.
 */

#include "SystemModelData.hpp"

#include <fstream>
#include <stdexcept>

#include <GridKit/Model/OptimalPowerFlow/SystemModelDataJSONParser.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    SystemModelData<double, size_t> parseSystemModelData(std::istream& stream)
    {
      return json::parse(stream).get<SystemModelData<double, size_t>>();
    }

    SystemModelData<double, size_t> parseSystemModelData(const std::filesystem::path& file)
    {
      std::ifstream stream(file);
      if (!stream)
      {
        throw std::runtime_error("Could not open file " + file.string());
      }
      return parseSystemModelData(stream);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
