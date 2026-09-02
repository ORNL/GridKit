#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser function implementation for the `VectorFitData` type
    ///
    /// The layout matches the fitting application's rational model files:
    /// `rows`, `cols`, `D` and `E` as 3x3 arrays, `poles` as an array of
    /// `[re, im]` pairs, and `residues` as `residues[q][row][col] = [re, im]`.
    /// See the `INPUT_FORMAT.md` in `GridKit/Model/EMT` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, VectorFitData<RealT, IdxT>& vf)
    {
      if (j.contains("rows"))
      {
        j.at("rows").get_to(vf.rows);
      }

      if (j.contains("cols"))
      {
        j.at("cols").get_to(vf.cols);
      }

      auto parse_matrix = [](const json& value, ABCMatrix<RealT>& matrix)
      {
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            value.at(n).at(k).get_to(matrix[n][k]);
          }
        }
      };

      if (j.contains("D"))
      {
        parse_matrix(j.at("D"), vf.D);
      }

      if (j.contains("E"))
      {
        parse_matrix(j.at("E"), vf.E);
      }

      if (j.contains("poles"))
      {
        for (const auto& raw_pole : j.at("poles"))
        {
          vf.poles.emplace_back(raw_pole.at(0).template get<RealT>(),
                                raw_pole.at(1).template get<RealT>());
        }
      }

      if (j.contains("residues"))
      {
        for (const auto& raw_residue : j.at("residues"))
        {
          auto& residue = vf.residues.emplace_back();
          for (size_t n = 0; n < 3; ++n)
          {
            for (size_t k = 0; k < 3; ++k)
            {
              residue[n][k] = {raw_residue.at(n).at(k).at(0).template get<RealT>(),
                               raw_residue.at(n).at(k).at(1).template get<RealT>()};
            }
          }
        }
      }
    }

    /**
     * @brief Parse a rational operator operand from a submodel JSON value.
     *
     * A string names a sidecar fit file, resolved against the working
     * directory unless absolute; an object is parsed inline.
     */
    template <typename RealT, typename IdxT>
    VectorFitData<RealT, IdxT> parseVectorFitOperand(const json& value)
    {
      if (value.is_string())
      {
        const std::filesystem::path file_path = value.template get<std::string>();
        std::ifstream               stream(file_path);
        if (!stream)
        {
          std::stringstream message;
          message << "Could not open fit file: " << file_path;
          Log::error() << message.str() << std::endl;
          throw std::runtime_error(message.str());
        }
        VectorFitData<RealT, IdxT> data(json::parse(stream));
        return data;
      }

      VectorFitData<RealT, IdxT> data(value);
      return data;
    }
  } // namespace EMT
} // namespace GridKit
