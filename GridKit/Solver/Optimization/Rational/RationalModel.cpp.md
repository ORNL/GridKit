```cpp
/**
 * @file RationalModel.cpp
 *
 * @brief Rational model evaluation and coefficient-contract validation.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "RationalModel.hpp"

#include <cmath>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief One residue matrix element for pole q.
     *
     * Residues are stored pole-major, row-major within each pole, so the
     * element (q, row, col) sits at (q * rows + row) * cols + col.
     */
    template <typename scalar_type, typename index_type>
    typename RationalModel<scalar_type, index_type>::ComplexT
    RationalModel<scalar_type, index_type>::residue(IdxT q,
                                                    IdxT row,
                                                    IdxT col) const
    {
      return residues[static_cast<size_t>((q * rows + row) * cols + col)];
    }

    /**
     * @brief Evaluate one matrix element of the model at j omega.
     *
     * Sums D + j omega E + sum over q of R_q / (j omega - p_q) directly;
     * correctness does not depend on conjugate ordering. An empty d or e
     * vector contributes zero. Evaluating exactly at a pole on the
     * imaginary axis yields a nonfinite value by design.
     */
    template <typename scalar_type, typename index_type>
    typename RationalModel<scalar_type, index_type>::ComplexT
    RationalModel<scalar_type, index_type>::evaluate(RealT omega,
                                                     IdxT  row,
                                                     IdxT  col) const
    {
      const auto     entry = static_cast<size_t>(row * cols + col);
      const ComplexT s{RealT{0}, omega};

      ComplexT value{RealT{0}, RealT{0}};
      if (!d.empty())
      {
        value += d[entry];
      }
      if (!e.empty())
      {
        value += s * e[entry];
      }

      const auto pole_count = static_cast<IdxT>(poles.size());
      for (IdxT q = 0; q < pole_count; ++q)
      {
        value += residue(q, row, col) / (s - poles[static_cast<size_t>(q)]);
      }
      return value;
    }

    /**
     * @brief Whether every pole satisfies Re(p) <= -margin.
     *
     * Computed from the poles; never a cached flag. An empty pole set is
     * stable.
     */
    template <typename scalar_type, typename index_type>
    bool RationalModel<scalar_type, index_type>::isStable(RealT margin) const
    {
      for (const auto& pole : poles)
      {
        if (pole.real() > -margin)
        {
          return false;
        }
      }
      return true;
    }

    /**
     * @brief Check the coefficient contract.
     *
     * Verifies storage sizes, finiteness of every coefficient, and the
     * conjugate-pair ordering: a pole is real when |Im p| is within the
     * scaled tolerance, and must then carry real residues; any other pole
     * must be followed immediately by its conjugate, residues likewise.
     * All comparisons are relative, |a - b| <= tolerance * (1 + |b|), so
     * zero tolerance demands exact conjugacy.
     *
     * @return 0 when valid, -1 on inconsistent sizes, -2 on a nonfinite
     *         value, -3 on a conjugate-pair violation
     */
    template <typename scalar_type, typename index_type>
    int RationalModel<scalar_type, index_type>::validate(RealT tolerance) const
    {
      if (rows < 1 || cols < 1)
      {
        return -1;
      }

      const auto channels =
          static_cast<size_t>(rows) * static_cast<size_t>(cols);
      if (residues.size() != poles.size() * channels)
      {
        return -1;
      }
      if (!d.empty() && d.size() != channels)
      {
        return -1;
      }
      if (!e.empty() && e.size() != channels)
      {
        return -1;
      }

      const auto finite = [](ComplexT value)
      {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
      };
      for (const auto& pole : poles)
      {
        if (!finite(pole))
        {
          return -2;
        }
      }
      for (const auto& value : residues)
      {
        if (!finite(value))
        {
          return -2;
        }
      }
      for (const auto value : d)
      {
        if (!std::isfinite(value))
        {
          return -2;
        }
      }
      for (const auto value : e)
      {
        if (!std::isfinite(value))
        {
          return -2;
        }
      }

      const auto real_pole = [tolerance](ComplexT pole)
      {
        return std::abs(pole.imag()) <= tolerance * (RealT{1} + std::abs(pole));
      };
      const auto conjugate = [tolerance](ComplexT value, ComplexT partner)
      {
        return std::abs(value - std::conj(partner)) <=
               tolerance * (RealT{1} + std::abs(partner));
      };

      const auto pole_count = poles.size();

      size_t q = 0;
      while (q < pole_count)
      {
        if (real_pole(poles[q]))
        {
          for (size_t entry = 0; entry < channels; ++entry)
          {
            const ComplexT value = residues[q * channels + entry];
            if (std::abs(value.imag()) >
                tolerance * (RealT{1} + std::abs(value)))
            {
              return -3;
            }
          }
          q += 1;
          continue;
        }

        if (q + 1 == pole_count)
        {
          return -3;
        }
        if (!conjugate(poles[q + 1], poles[q]))
        {
          return -3;
        }
        for (size_t entry = 0; entry < channels; ++entry)
        {
          if (!conjugate(residues[(q + 1) * channels + entry],
                         residues[q * channels + entry]))
          {
            return -3;
          }
        }
        q += 2;
      }

      return 0;
    }

    template struct RationalModel<double, long int>;
    template struct RationalModel<double, size_t>;
  } // namespace Optimization
} // namespace GridKit
```
