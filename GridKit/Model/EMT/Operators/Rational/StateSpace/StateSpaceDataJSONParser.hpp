#pragma once

#include <filesystem>
#include <fstream>

#include <GridKit/Model/EMT/Operators/Rational/RationalDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/StateSpace/StateSpaceData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename RealT, typename IdxT>
    void from_json(const nlohmann::json& j, StateSpaceData<RealT, IdxT>& result)
    {
      StateSpaceData<RealT, IdxT> data;
      data.rows = parseRationalDimension<IdxT>(j, "rows");
      data.cols = parseRationalDimension<IdxT>(j, "cols");
      if (data.rows <= 0 || data.cols <= 0)
      {
        throw std::invalid_argument("StateSpace: dimensions must be positive");
      }
      const size_t rows = static_cast<size_t>(data.rows), cols = static_cast<size_t>(data.cols);
      data.D    = RationalMatrix<RealT>(rows, cols);
      data.E    = RationalMatrix<RealT>(rows, cols);
      auto real = [](const nlohmann::json& value)
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
          throw std::invalid_argument("StateSpace: poles must be an array");
        }
        for (const auto& pole : j.at("poles"))
        {
          data.poles.push_back(parseRationalComplex<RealT>(pole));
        }
      }
      const size_t Q = data.poles.size();
      data.C         = RationalMatrix<std::complex<RealT>>(rows, Q);
      data.B         = RationalMatrix<std::complex<RealT>>(Q, cols);
      if (Q && (!j.contains("B") || !j.contains("C")))
      {
        throw std::invalid_argument("StateSpace: poles require B and C factors");
      }
      if (j.contains("C"))
      {
        data.C = parseRationalMatrix<std::complex<RealT>>(j.at("C"), rows, Q, parseRationalComplex<RealT>);
      }
      if (j.contains("B"))
      {
        data.B = parseRationalMatrix<std::complex<RealT>>(j.at("B"), Q, cols, parseRationalComplex<RealT>);
      }
      if (data.validate())
      {
        throw std::invalid_argument("StateSpace: invalid rational coefficients");
      }
      result = std::move(data);
    }

    template <typename RealT, typename IdxT>
    StateSpaceData<RealT, IdxT> parseStateSpaceOperand(const nlohmann::json& value)
    {
      if (value.is_string())
      {
        const std::filesystem::path path = value.template get<std::string>();
        std::ifstream               stream(path);
        if (!stream)
        {
          throw std::runtime_error("Could not open fit file: " + path.string());
        }
        return nlohmann::json::parse(stream).template get<StateSpaceData<RealT, IdxT>>();
      }
      return value.template get<StateSpaceData<RealT, IdxT>>();
    }
  } // namespace EMT
} // namespace GridKit
