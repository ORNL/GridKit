#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Operators/Rational/RationalDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// Parse matrix dimensions and pole-major residue data without truncation.
    template <typename RealT, typename IdxT>
    void from_json(const json& j, VectorFitData<RealT, IdxT>& vf)
    {
      VectorFitData<RealT, IdxT> data;
      data.rows = parseRationalDimension<IdxT>(j, "rows");
      data.cols = parseRationalDimension<IdxT>(j, "cols");
      if (data.rows <= 0 || data.cols <= 0)
      {
        throw std::invalid_argument("VectorFit: dimensions must be positive");
      }
      const size_t rows = static_cast<size_t>(data.rows), cols = static_cast<size_t>(data.cols);
      data.D    = RationalMatrix<RealT>(rows, cols);
      data.E    = RationalMatrix<RealT>(rows, cols);
      auto real = [](const json& value)
      { return value.template get<RealT>(); };
      if (j.contains("D"))
      {
        data.D = parseRationalMatrix<RealT>(j.at("D"), rows, cols, real);
      }
      if (j.contains("E"))
      {
        data.E = parseRationalMatrix<RealT>(j.at("E"), rows, cols, real);
      }
      if (j.contains("poles"))
      {
        if (!j.at("poles").is_array())
        {
          throw std::invalid_argument("VectorFit: poles must be an array");
        }
        for (const auto& pole : j.at("poles"))
        {
          data.poles.push_back(parseRationalComplex<RealT>(pole));
        }
      }
      if (j.contains("residues"))
      {
        if (!j.at("residues").is_array())
        {
          throw std::invalid_argument("VectorFit: residues must be an array");
        }
        for (const auto& residue : j.at("residues"))
        {
          data.residues.push_back(parseRationalMatrix<std::complex<RealT>>(residue, rows, cols, parseRationalComplex<RealT>));
        }
      }
      if (data.validate())
      {
        throw std::invalid_argument("VectorFit: invalid rational coefficients");
      }
      vf = std::move(data);
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
