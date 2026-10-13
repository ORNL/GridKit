#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace GridKit::Examples::EMT
{
  using Complex    = std::complex<double>;
  using Expression = std::vector<std::pair<size_t, double>>;

  /// Sparse linear equations K*y + M*y' = 0. Assembly occurs only at startup.
  struct LinearDAE
  {
    struct Entry
    {
      size_t row, col;
      double k, m;
    };

    std::vector<bool>   differential;
    std::vector<double> scale;
    std::vector<Entry>  entries;

    size_t variables(size_t count, bool diff, double magnitude)
    {
      const size_t first = scale.size();
      differential.insert(differential.end(), count, diff);
      scale.insert(scale.end(), count, magnitude);
      return first;
    }

    void add(size_t row, size_t col, double k, double m = 0)
    {
      entries.push_back({row, col, k, m});
    }

    void add(size_t row, const Expression& expression, double factor)
    {
      for (const auto& [col, value] : expression)
      {
        add(row, col, factor * value);
      }
    }

    void finish()
    {
      std::map<std::pair<size_t, size_t>, std::pair<double, double>> summed;
      for (const auto& e : entries)
      {
        auto& [k, m]  = summed[{e.row, e.col}];
        k            += e.k;
        m            += e.m;
      }
      entries.clear();
      for (const auto& [index, value] : summed)
      {
        entries.push_back({index.first, index.second, value.first, value.second});
      }
    }
  };

  inline double evaluate(const Expression& expression, std::span<const double> values)
  {
    double result = 0;
    for (const auto& [index, coefficient] : expression)
    {
      result += coefficient * values[index];
    }
    return result;
  }

  /// Small dense complex solve used only for sinusoidal initialization.
  inline std::vector<Complex> solve(std::vector<Complex> matrix, std::vector<Complex> rhs)
  {
    const size_t n = rhs.size();
    for (size_t row = 0; row < n; ++row)
    {
      double scale = 0;
      for (size_t col = 0; col < n; ++col)
      {
        scale = std::max(scale, std::abs(matrix[row * n + col]));
      }
      if (!(scale > 0))
      {
        throw std::runtime_error("Singular initialization row");
      }
      for (size_t col = 0; col < n; ++col)
      {
        matrix[row * n + col] /= scale;
      }
      rhs[row] /= scale;
    }
    for (size_t col = 0; col < n; ++col)
    {
      size_t pivot = col;
      for (size_t row = col + 1; row < n; ++row)
      {
        if (std::abs(matrix[row * n + col]) > std::abs(matrix[pivot * n + col]))
        {
          pivot = row;
        }
      }
      if (std::abs(matrix[pivot * n + col]) < 1e-14)
      {
        throw std::runtime_error("Singular sinusoidal circuit");
      }
      for (size_t j = 0; j < n; ++j)
      {
        std::swap(matrix[col * n + j], matrix[pivot * n + j]);
      }
      std::swap(rhs[col], rhs[pivot]);
      for (size_t row = col + 1; row < n; ++row)
      {
        const Complex factor = matrix[row * n + col] / matrix[col * n + col];
        for (size_t j = col + 1; j < n; ++j)
        {
          matrix[row * n + j] -= factor * matrix[col * n + j];
        }
        rhs[row] -= factor * rhs[col];
      }
    }
    for (size_t row = n; row-- > 0;)
    {
      for (size_t col = row + 1; col < n; ++col)
      {
        rhs[row] -= matrix[row * n + col] * rhs[col];
      }
      rhs[row] /= matrix[row * n + row];
    }
    return rhs;
  }
} // namespace GridKit::Examples::EMT
