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
    /// Factorized coefficients H(s) = D + s E + C (s I - P)^-1 B.
    template <typename real_type, typename index_type>
    struct StateSpaceData
    {
      using RealT    = real_type;
      using IdxT     = index_type;
      using ComplexT = std::complex<RealT>;
      IdxT                     rows{3}, cols{3};
      RationalMatrix<RealT>    D, E;
      std::vector<ComplexT>    poles;
      RationalMatrix<ComplexT> C{3, 0}, B{0, 3};

      int validate() const
      {
        if (rows <= 0 || cols <= 0)
        {
          return 1;
        }
        const size_t N = static_cast<size_t>(rows), K = static_cast<size_t>(cols), Q = poles.size();
        int          errors = 0;
        if (!D.hasShape(N, K) || !E.hasShape(N, K) || !C.hasShape(N, Q) || !B.hasShape(Q, K))
        {
          return 1;
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
        for (const auto* matrix : {&C, &B})
        {
          for (const auto& row : *matrix)
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
        for (size_t q = 0; q < Q; ++q)
        {
          if (!std::isfinite(poles[q].real()) || !std::isfinite(poles[q].imag()))
          {
            ++errors;
            continue;
          }
          if (poles[q].real() >= RealT{0})
          {
            ::GridKit::Utilities::Logger::warning() << "StateSpace: pole " << q << " is not stable\n";
          }
          if (poles[q].imag() == RealT{0})
          {
            for (size_t n = 0; n < N; ++n)
            {
              if (C[n][q].imag() != RealT{0})
              {
                ++errors;
              }
            }
            for (size_t k = 0; k < K; ++k)
            {
              if (B[q][k].imag() != RealT{0})
              {
                ++errors;
              }
            }
          }
          else if (q + 1 >= Q || poles[q + 1] != std::conj(poles[q]))
          {
            ++errors;
          }
          else
          {
            for (size_t n = 0; n < N; ++n)
            {
              if (C[n][q + 1] != std::conj(C[n][q]))
              {
                ++errors;
              }
            }
            for (size_t k = 0; k < K; ++k)
            {
              if (B[q + 1][k] != std::conj(B[q][k]))
              {
                ++errors;
              }
            }
            ++q;
          }
        }
        return errors;
      }
    };
  } // namespace EMT
} // namespace GridKit
