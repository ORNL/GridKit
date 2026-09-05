#pragma once

#include <complex>
#include <limits>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Operators/Rational/RationalMatrix.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename IdxT>
    IdxT parseRationalDimension(const nlohmann::json& object, const char* key)
    {
      if (!object.contains(key))
      {
        return IdxT{3};
      }
      const auto& value = object.at(key);
      if (!value.is_number_integer() || (value.is_number_integer() && !value.is_number_unsigned() && value.template get<int64_t>() <= 0))
      {
        throw std::invalid_argument("Rational: dimensions must be positive integers");
      }
      const auto dimension = value.template get<uint64_t>();
      if (dimension == 0 || dimension > static_cast<uint64_t>(std::numeric_limits<IdxT>::max()))
      {
        throw std::invalid_argument("Rational: dimension exceeds the index range");
      }
      return static_cast<IdxT>(dimension);
    }

    template <typename RealT>
    std::complex<RealT> parseRationalComplex(const nlohmann::json& value)
    {
      if (!value.is_array() || value.size() != 2)
      {
        throw std::invalid_argument("Rational: complex values require [real, imaginary]");
      }
      return {value.at(0).template get<RealT>(), value.at(1).template get<RealT>()};
    }

    template <typename T, typename Parse>
    RationalMatrix<T> parseRationalMatrix(const nlohmann::json& value, size_t rows, size_t cols, Parse parse)
    {
      if (!value.is_array() || value.size() != rows)
      {
        throw std::invalid_argument("Rational: matrix row count mismatch");
      }
      RationalMatrix<T> matrix(rows, cols);
      for (size_t n = 0; n < rows; ++n)
      {
        if (!value.at(n).is_array() || value.at(n).size() != cols)
        {
          throw std::invalid_argument("Rational: matrix column count mismatch");
        }
        for (size_t k = 0; k < cols; ++k)
        {
          matrix[n][k] = parse(value.at(n).at(k));
        }
      }
      return matrix;
    }
  } // namespace EMT
} // namespace GridKit
