#pragma once

#include <filesystem>
#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

#include "LinearDAE.hpp"

namespace GridKit::Examples::EMT
{
  using Json = nlohmann::json;

  inline void fields(const Json& j, std::initializer_list<const char*> allowed)
  {
    if (!j.is_object())
    {
      throw std::invalid_argument("Expected a JSON object");
    }
    for (const auto& [key, value] : j.items())
    {
      (void) value;
      if (std::none_of(allowed.begin(), allowed.end(), [&](const char* name)
                       { return key == name; }))
      {
        throw std::invalid_argument("Unknown field: " + key);
      }
    }
  }

  inline double real(const Json& j)
  {
    if (!j.is_number() || !std::isfinite(j.get<double>()))
    {
      throw std::invalid_argument("Expected a finite real number");
    }
    return j.get<double>();
  }

  inline size_t dimension(const Json& j)
  {
    if (!j.is_number_integer() || j.get<long long>() <= 0)
    {
      throw std::invalid_argument("Expected a positive integer dimension");
    }
    return j.get<size_t>();
  }

  inline Complex complex(const Json& j)
  {
    if (j.is_number())
    {
      return {real(j), 0};
    }
    if (!j.is_array() || j.size() != 2)
    {
      throw std::invalid_argument("Complex coefficient must be [real, imaginary]");
    }
    return {real(j[0]), real(j[1])};
  }

  inline Json readJson(const std::filesystem::path& path)
  {
    std::ifstream stream(path);
    if (!stream)
    {
      throw std::runtime_error("Cannot open " + path.string());
    }
    return Json::parse(stream);
  }

  struct VectorFit
  {
    size_t                            rows = 0, cols = 0;
    std::vector<double>               D;
    std::vector<Complex>              poles;
    std::vector<std::vector<Complex>> residues;

    static VectorFit identity(size_t n)
    {
      VectorFit fit;
      fit.rows = fit.cols = n;
      fit.D.assign(n * n, 0);
      for (size_t i = 0; i < n; ++i)
      {
        fit.D[i * n + i] = 1;
      }
      return fit;
    }

    static VectorFit parse(const Json& value, const std::filesystem::path& directory)
    {
      if (value.is_string())
      {
        const auto path = directory / value.get<std::string>();
        return parse(readJson(path), path.parent_path());
      }
      fields(value, {"rows", "cols", "D", "E", "poles", "residues"});
      VectorFit fit;
      fit.rows    = dimension(value.at("rows"));
      fit.cols    = dimension(value.at("cols"));
      auto matrix = [&](const Json& j)
      {
        std::vector<Complex> result;
        if (!j.is_array() || j.size() != fit.rows)
        {
          throw std::invalid_argument("VectorFit matrix row count");
        }
        for (const auto& row : j)
        {
          if (!row.is_array() || row.size() != fit.cols)
          {
            throw std::invalid_argument("VectorFit matrix column count");
          }
          for (const auto& item : row)
          {
            result.push_back(complex(item));
          }
        }
        return result;
      };
      fit.D.assign(fit.rows * fit.cols, 0);
      if (value.contains("D"))
      {
        const auto d = matrix(value.at("D"));
        for (size_t i = 0; i < d.size(); ++i)
        {
          if (d[i].imag() != 0)
          {
            throw std::invalid_argument("VectorFit D must be real");
          }
          fit.D[i] = d[i].real();
        }
      }
      if (value.contains("E"))
      {
        for (Complex e : matrix(value.at("E")))
        {
          if (e != Complex{})
          {
            throw std::invalid_argument("LineDistributed requires proper fits: E = 0");
          }
        }
      }
      const auto p = value.value("poles", Json::array());
      const auto r = value.value("residues", Json::array());
      if (!p.is_array() || !r.is_array() || p.size() != r.size())
      {
        throw std::invalid_argument("VectorFit pole/residue count");
      }
      for (const auto& item : p)
      {
        fit.poles.push_back(complex(item));
        if (!(fit.poles.back().real() < 0))
        {
          throw std::invalid_argument("VectorFit requires strictly stable poles");
        }
      }
      for (const auto& item : r)
      {
        fit.residues.push_back(matrix(item));
      }
      auto conjugate = [](Complex a, Complex b)
      { return std::abs(a - std::conj(b)) <= 1e-12 * std::max({1.0, std::abs(a), std::abs(b)}); };
      for (size_t q = 0; q < p.size(); ++q)
      {
        if (fit.poles[q].imag() == 0)
        {
          for (Complex residue : fit.residues[q])
          {
            if (residue.imag() != 0)
            {
              throw std::invalid_argument("Real pole requires real residues");
            }
          }
        }
        else
        {
          if (q + 1 == p.size() || !conjugate(fit.poles[q], fit.poles[q + 1]))
          {
            throw std::invalid_argument("Poles must have adjacent conjugate partners");
          }
          for (size_t i = 0; i < fit.D.size(); ++i)
          {
            if (!conjugate(fit.residues[q][i], fit.residues[q + 1][i]))
            {
              throw std::invalid_argument("Residues must have conjugate partners");
            }
          }
          ++q;
        }
      }
      return fit;
    }

    size_t allocate(LinearDAE& dae, double input_scale) const
    {
      return dae.variables(cols * poles.size(), true, input_scale);
    }

    /// Real pole blocks. x = |p|*w normalizes memory states to input units.
    /// Complex pairs occupy 2*K real states and use both residue matrices.
    std::vector<Expression> stamp(LinearDAE& dae, size_t state, size_t input, std::span<const size_t> mapping = {}) const
    {
      auto column = [&](size_t k)
      { return input + (mapping.empty() ? k : mapping[k]); };
      std::vector<Expression> output(rows);
      for (size_t n = 0; n < rows; ++n)
      {
        for (size_t k = 0; k < cols; ++k)
        {
          if (D[n * cols + k] != 0)
          {
            output[n].push_back({column(k), D[n * cols + k]});
          }
        }
      }
      for (size_t q = 0; q < poles.size(); ++q)
      {
        const double a = poles[q].real(), b = poles[q].imag();
        const double scale = std::abs(poles[q]);
        const bool   pair  = b != 0;
        const size_t w = state + q * cols, v = w + cols;
        for (size_t k = 0; k < cols; ++k)
        {
          dae.add(w + k, w + k, a, -1);
          dae.add(w + k, column(k), scale);
          if (pair)
          {
            dae.add(w + k, v + k, -b);
            dae.add(v + k, w + k, b);
            dae.add(v + k, v + k, a, -1);
          }
        }
        for (size_t n = 0; n < rows; ++n)
        {
          for (size_t k = 0; k < cols; ++k)
          {
            const Complex c = (pair ? 2.0 : 1.0) * residues[q][n * cols + k] / scale;
            if (c.real() != 0)
            {
              output[n].push_back({w + k, c.real()});
            }
            if (pair && c.imag() != 0)
            {
              output[n].push_back({v + k, -c.imag()});
            }
          }
        }
        if (pair)
        {
          ++q;
        }
      }
      return output;
    }

    std::vector<Complex> frequency(Complex s) const
    {
      std::vector<Complex> result(D.begin(), D.end());
      for (size_t q = 0; q < poles.size(); ++q)
      {
        for (size_t i = 0; i < result.size(); ++i)
        {
          result[i] += residues[q][i] / (s - poles[q]);
        }
      }
      return result;
    }
  };
} // namespace GridKit::Examples::EMT
