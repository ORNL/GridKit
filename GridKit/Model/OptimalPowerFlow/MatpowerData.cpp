/**
 * @file MatpowerData.cpp
 * @brief MATPOWER case file input.
 */

#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>

#include <GridKit/Model/OptimalPowerFlow/MatpowerData.hpp>
#include <GridKit/Model/PowerFlow/MatpowerParser.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    const MatpowerMatrix& MatpowerData::matrix(const std::string& name) const
    {
      const auto entry = matrices.find(name);
      if (entry == matrices.end())
      {
        throw std::invalid_argument("MATPOWER case has no mpc." + name);
      }
      return entry->second;
    }

    MatpowerData parseMatpowerData(std::istream& stream)
    {
      static const std::regex matrix_start(R"(\s*mpc\.(\w+)\s*=\s*\[(.*))");

      MatpowerData data;
      for (std::string line; std::getline(stream, line);)
      {
        std::smatch match;
        if (!std::regex_match(line, match, matrix_start))
        {
          continue;
        }

        auto& matrix = data.matrices[match[1]];
        GridKit::readMatPowerMatrix(stream, line, [&](const std::string& row)
                                    {
          std::istringstream  values(row);
          std::vector<double> numbers{std::istream_iterator<double>(values), std::istream_iterator<double>()};
          if (!values.eof())
          {
            throw std::invalid_argument("MATPOWER matrix row is not numeric: " + row);
          }
          matrix.push_back(std::move(numbers)); });
      }
      return data;
    }

    MatpowerData parseMatpowerData(const std::filesystem::path& file)
    {
      std::ifstream stream(file);
      if (!stream)
      {
        throw std::runtime_error("Could not open MATPOWER case file " + file.string());
      }
      return parseMatpowerData(stream);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
