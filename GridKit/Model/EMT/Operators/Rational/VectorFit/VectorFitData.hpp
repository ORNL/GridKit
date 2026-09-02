/**
 * @file VectorFitData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT vector-fitted rational operators
 *
 */
#pragma once

#include <cmath>
#include <complex>
#include <vector>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Contains modeling data for a vector-fitted rational operator
     *
     * The operator realizes H(s) = D + sE + sum_q R_q / (s - p_q). Complex
     * values appear in data only; the runtime model splits them into real
     * derived parameters at construction.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct VectorFitData
    {
      using RealT    = real_type;
      using IdxT     = index_type;
      using ComplexT = std::complex<RealT>;

      IdxT rows{3}; ///< Output dimension
      IdxT cols{3}; ///< Input dimension

      ABCMatrix<RealT> D{}; ///< Constant coefficient
      ABCMatrix<RealT> E{}; ///< Linear coefficient

      std::vector<ComplexT>            poles;    ///< Poles, conjugate pairs adjacent
      std::vector<ABCMatrix<ComplexT>> residues; ///< Residue matrices, pole-major

      /**
       * @brief Check the data against the model specification.
       *
       * Real poles carry real residues, and each nonreal pole and residue is
       * followed by its exact conjugate. An unstable pole is reported as a
       * warning, not an error.
       *
       * @return Number of specification violations found
       */
      int validate() const
      {
        using Log = ::GridKit::Utilities::Logger;

        int error_count = 0;

        if (rows != 3 || cols != 3)
        {
          Log::error() << "VectorFit: the input and output dimensions must be three\n";
          ++error_count;
        }

        if (residues.size() != poles.size())
        {
          Log::error() << "VectorFit: each pole requires one residue matrix\n";
          ++error_count;
          return error_count;
        }

        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            if (!std::isfinite(D[n][k]) || !std::isfinite(E[n][k]))
            {
              Log::error() << "VectorFit: the coefficient matrices must be finite\n";
              ++error_count;
            }
          }
        }

        size_t q = 0;
        while (q < poles.size())
        {
          if (!std::isfinite(poles[q].real()) || !std::isfinite(poles[q].imag()))
          {
            Log::error() << "VectorFit: the poles must be finite\n";
            ++error_count;
            ++q;
            continue;
          }

          if (poles[q].real() >= 0.0)
          {
            Log::warning() << "VectorFit: pole " << q << " is not stable\n";
          }

          if (poles[q].imag() == 0.0)
          {
            for (size_t n = 0; n < 3; ++n)
            {
              for (size_t k = 0; k < 3; ++k)
              {
                if (residues[q][n][k].imag() != 0.0)
                {
                  Log::error() << "VectorFit: the residue of real pole " << q
                               << " must be real\n";
                  ++error_count;
                }
              }
            }
            ++q;
            continue;
          }

          if (q + 1 >= poles.size() || poles[q + 1] != std::conj(poles[q]))
          {
            Log::error() << "VectorFit: nonreal pole " << q
                         << " must be followed by its conjugate\n";
            ++error_count;
            ++q;
            continue;
          }

          for (size_t n = 0; n < 3; ++n)
          {
            for (size_t k = 0; k < 3; ++k)
            {
              if (residues[q + 1][n][k] != std::conj(residues[q][n][k]))
              {
                Log::error() << "VectorFit: the residue of pole " << q + 1
                             << " must conjugate the residue of pole " << q << "\n";
                ++error_count;
              }
            }
          }
          q += 2;
        }

        return error_count;
      }
    };
  } // namespace EMT
} // namespace GridKit
