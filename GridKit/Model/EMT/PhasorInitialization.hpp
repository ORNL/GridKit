#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit::EMT
{
  /** Whether the nonzero derivative columns define independent current variables. */
  template <typename RealT, typename MatrixT>
  bool independentDerivativeColumns(const MatrixT& coefficient)
  {
    ABCMatrix<RealT> matrix{};
    RealT            scale   = 0;
    size_t           columns = 0;
    for (size_t k = 0; k < 3; ++k)
    {
      bool nonzero = false;
      for (size_t n = 0; n < 3; ++n)
      {
        if (!std::isfinite(coefficient[n][k]))
          return false;
        nonzero = nonzero || coefficient[n][k] != RealT{0};
        scale   = std::max(scale, std::abs(coefficient[n][k]));
      }
      if (nonzero)
      {
        for (size_t n = 0; n < 3; ++n)
          matrix[n][columns] = coefficient[n][k];
        ++columns;
      }
    }
    if (scale != RealT{0})
    {
      for (auto& row : matrix)
      {
        for (auto& value : row)
          value /= scale;
      }
    }
    const RealT tolerance = 3 * std::numeric_limits<RealT>::epsilon();
    for (size_t k = 0; k < columns; ++k)
    {
      size_t pivot = k;
      for (size_t n = k + 1; n < 3; ++n)
        if (std::abs(matrix[n][k]) > std::abs(matrix[pivot][k]))
          pivot = n;
      if (std::abs(matrix[pivot][k]) <= tolerance)
        return false;
      std::swap(matrix[k], matrix[pivot]);
      for (size_t n = k + 1; n < 3; ++n)
      {
        const RealT multiplier = matrix[n][k] / matrix[k][k];
        for (size_t j = k + 1; j < columns; ++j)
          matrix[n][j] -= multiplier * matrix[k][j];
      }
    }
    return true;
  }

  /**
   * @brief Solve H(j omega) i = v from instantaneous sinusoidal samples.
   *
   * Partial pivoting supports coupled phase matrices. A singular transfer or
   * nonfinite input returns an error without changing the output samples.
   */
  template <typename RealT>
  int solvePhasorSystem(const ABCMatrix<RealT>& H_re,
                        const ABCMatrix<RealT>& H_im,
                        RealT                   omega,
                        const ABCVector<RealT>& v,
                        const ABCVector<RealT>& v_dot,
                        ABCVector<RealT>&       i,
                        ABCVector<RealT>&       i_dot)
  {
    using ComplexT = std::complex<RealT>;
    if (!std::isfinite(omega) || omega <= RealT{0})
    {
      return 1;
    }

    ABCMatrix<ComplexT> matrix{};
    ABCVector<ComplexT> rhs{};
    for (size_t n = 0; n < 3; ++n)
    {
      if (!std::isfinite(v[n]) || !std::isfinite(v_dot[n]))
      {
        return 1;
      }
      rhs[n] = {v[n], -v_dot[n] / omega};
      for (size_t k = 0; k < 3; ++k)
      {
        if (!std::isfinite(H_re[n][k]) || !std::isfinite(H_im[n][k]))
        {
          return 1;
        }
        matrix[n][k] = {H_re[n][k], H_im[n][k]};
      }
    }

    for (size_t k = 0; k < 3; ++k)
    {
      size_t pivot = k;
      for (size_t n = k + 1; n < 3; ++n)
      {
        if (std::abs(matrix[n][k]) > std::abs(matrix[pivot][k]))
        {
          pivot = n;
        }
      }
      if (std::abs(matrix[pivot][k]) == RealT{0})
      {
        return 1;
      }
      std::swap(matrix[k], matrix[pivot]);
      std::swap(rhs[k], rhs[pivot]);
      for (size_t n = k + 1; n < 3; ++n)
      {
        const ComplexT multiplier = matrix[n][k] / matrix[k][k];
        for (size_t j = k + 1; j < 3; ++j)
        {
          matrix[n][j] -= multiplier * matrix[k][j];
        }
        rhs[n] -= multiplier * rhs[k];
      }
    }

    for (size_t n = 3; n-- > 0;)
    {
      for (size_t k = n + 1; k < 3; ++k)
      {
        rhs[n] -= matrix[n][k] * rhs[k];
      }
      rhs[n] /= matrix[n][n];
      if (!std::isfinite(rhs[n].real()) || !std::isfinite(rhs[n].imag())
          || !std::isfinite(-omega * rhs[n].imag()))
      {
        return 1;
      }
    }
    for (size_t n = 0; n < 3; ++n)
    {
      i[n]     = rhs[n].real();
      i_dot[n] = -omega * rhs[n].imag();
    }
    return 0;
  }
} // namespace GridKit::EMT
