#pragma once

#include <cmath>
#include <complex>
#include <vector>

#include <GridKit/Model/EMT/Operators/Rational/RationalMatrix.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Coefficients of H(s) = D + s E + sum R_q / (s - p_q).
    template <typename real_type, typename index_type>
    struct VectorFitData
    {
      using RealT    = real_type;
      using IdxT     = index_type;
      using ComplexT = std::complex<RealT>;
      IdxT                                  rows{3}, cols{3};
      RationalMatrix<RealT>                 D, E;
      std::vector<ComplexT>                 poles;
      std::vector<RationalMatrix<ComplexT>> residues;

      int validate() const
      {
        using Log  = ::GridKit::Utilities::Logger;
        int errors = 0;
        if (rows <= 0 || cols <= 0)
        {
          Log::error() << "VectorFit: input and output dimensions must be positive\n";
          return 1;
        }
        const size_t N = static_cast<size_t>(rows), K = static_cast<size_t>(cols);
        if (!D.hasShape(N, K) || !E.hasShape(N, K))
        {
          Log::error() << "VectorFit: coefficient dimensions do not match rows and cols\n";
          ++errors;
        }
        for (const auto* matrix : {&D, &E})
        {
          for (const auto& row : *matrix)
          {
            for (const auto value : row)
            {
              if (!std::isfinite(value))
              {
                ++errors;
              }
            }
          }
        }
        if (residues.size() != poles.size())
        {
          Log::error() << "VectorFit: each pole requires one residue matrix\n";
          return errors + 1;
        }
        for (const auto& matrix : residues)
        {
          if (!matrix.hasShape(N, K))
          {
            ++errors;
          }
          for (const auto& row : matrix)
          {
            for (const auto value : row)
            {
              if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
              {
                ++errors;
              }
            }
          }
        }
        if (errors)
        {
          return errors;
        }
        for (size_t q = 0; q < poles.size(); ++q)
        {
          const auto p = poles[q];
          if (!std::isfinite(p.real()) || !std::isfinite(p.imag()))
          {
            ++errors;
            continue;
          }
          if (p.real() >= RealT{0})
          {
            Log::warning() << "VectorFit: pole " << q << " is not stable\n";
          }
          if (p.imag() == RealT{0})
          {
            for (const auto& row : residues[q])
            {
              for (const auto value : row)
              {
                if (value.imag() != RealT{0})
                {
                  ++errors;
                }
              }
            }
          }
          else if (q + 1 >= poles.size() || poles[q + 1] != std::conj(p))
          {
            ++errors;
          }
          else
          {
            for (size_t n = 0; n < N; ++n)
            {
              for (size_t k = 0; k < K; ++k)
              {
                if (residues[q + 1][n][k] != std::conj(residues[q][n][k]))
                {
                  ++errors;
                }
              }
            }
            ++q;
          }
        }
        if (errors)
        {
          Log::error() << "VectorFit: invalid pole or residue data\n";
        }
        return errors;
      }
    };
  } // namespace EMT
} // namespace GridKit
