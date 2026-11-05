/**
 * @file SystemModelData.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT system data from a JSON case.
 */

#include "SystemModelData.hpp"

#include <fstream>
#include <stdexcept>

#include <GridKit/Model/EMT/SystemModelDataJSONParser.hpp>

namespace GridKit
{
  namespace EMT
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
  } // namespace EMT
} // namespace GridKit
