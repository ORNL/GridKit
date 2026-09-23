/**
 * @file MatpowerData.cpp
 * @brief MATPOWER case file input.
 */

#include <algorithm>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>

#include <GridKit/Model/OptimalPowerFlow/MatpowerData.hpp>

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

      MatpowerData    data;
      MatpowerMatrix* matrix = nullptr;
      for (std::string line; std::getline(stream, line);)
      {
        std::string text = line.substr(0, line.find('%'));

        std::smatch match;
        if (matrix == nullptr)
        {
          if (!std::regex_match(text, match, matrix_start))
          {
            continue;
          }
          matrix = &data.matrices[match[1]];
          text   = match[2];
        }

        // Rows end at `;` or at the end of the line, and `]` ends the matrix
        const size_t end = text.find(']');
        text             = text.substr(0, end);
        std::replace(text.begin(), text.end(), ';', '\n');

        std::istringstream rows(text);
        for (std::string row; std::getline(rows, row);)
        {
          std::istringstream  values(row);
          std::vector<double> numbers{std::istream_iterator<double>(values), std::istream_iterator<double>()};
          if (!values.eof())
          {
            throw std::invalid_argument("MATPOWER matrix row is not numeric: " + row);
          }
          if (!numbers.empty())
          {
            matrix->push_back(std::move(numbers));
          }
        }

        if (end != std::string::npos)
        {
          matrix = nullptr;
        }
      }

      if (matrix != nullptr)
      {
        throw std::invalid_argument("MATPOWER matrix is not closed by ]");
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
