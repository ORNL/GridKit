#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace GridKit::PhasorDynamics::detail
{
  template <typename RealT, size_t N>
  struct HistoryPolynomial
  {
    RealT                begin;
    RealT                end;
    std::array<RealT, N> coefficients;

    RealT evaluate(RealT t) const
    {
      auto        c = coefficients;
      const RealT x = (t - begin) / (end - begin);
      for (size_t n = N - 1; n > 0; --n)
      {
        for (size_t i = 0; i < n; ++i)
        {
          c[i] = (1 - x) * c[i] + x * c[i + 1];
        }
      }
      return c[0];
    }

    std::pair<HistoryPolynomial, HistoryPolynomial> split(RealT t) const
    {
      auto              c = coefficients;
      HistoryPolynomial left{begin, t, {}}, right{t, end, {}};
      const RealT       x       = (t - begin) / (end - begin);
      left.coefficients[0]      = c[0];
      right.coefficients[N - 1] = c[N - 1];
      for (size_t n = N - 1; n > 0; --n)
      {
        for (size_t i = 0; i < n; ++i)
        {
          c[i] = (1 - x) * c[i] + x * c[i + 1];
        }
        left.coefficients[N - n]  = c[0];
        right.coefficients[n - 1] = c[n - 1];
      }
      return {left, right};
    }

    HistoryPolynomial restrictTo(RealT a, RealT b) const
    {
      auto result = *this;
      if (b < result.end)
      {
        result = result.split(b).first;
      }
      if (a > result.begin)
      {
        result = result.split(a).second;
      }
      return result;
    }

    /// Reconstruct a quintic at the six Chebyshev-Lobatto nodes.
    static HistoryPolynomial interpolate(RealT begin, RealT end, const std::array<RealT, N>& values)
    {
      static_assert(N == 6);
      // Inverse Bernstein collocation matrix: V_ij = binomial(5,j) x_i^j (1-x_i)^(5-j).
      static constexpr RealT inverse[6][6] = {
          {1, 0, 0, 0, 0, 0},
          {-2.4, 4.1888543819998318, -1.1577708763999663, 0.6111456180001682, -0.44222912360003364, 0.2},
          {2.52, -5.2944271909999161, 5.8499689437998486, -3.5055728090000842, 2.6300310562001514, -1.2},
          {-1.2, 2.6300310562001514, -3.5055728090000842, 5.8499689437998486, -5.2944271909999161, 2.52},
          {0.2, -0.44222912360003364, 0.6111456180001682, -1.1577708763999663, 4.1888543819998318, -2.4},
          {0, 0, 0, 0, 0, 1}};
      HistoryPolynomial output{begin, end, {}};
      for (size_t i = 0; i < 6; ++i)
      {
        for (size_t j = 0; j < 6; ++j)
        {
          output.coefficients[i] += inverse[i][j] * values[j];
        }
      }
      return output;
    }

    template <size_t M>
    HistoryPolynomial<RealT, M> elevate() const
    {
      static_assert(M >= N);
      HistoryPolynomial<RealT, M> result{begin, end, {}};
      std::copy(coefficients.begin(), coefficients.end(), result.coefficients.begin());
      for (size_t degree = N; degree < M; ++degree)
      {
        result.coefficients[degree] = result.coefficients[degree - 1];
        for (size_t i = degree - 1; i > 0; --i)
        {
          const RealT weight     = static_cast<RealT>(i) / static_cast<RealT>(degree);
          result.coefficients[i] = weight * result.coefficients[i - 1]
                                   + (1 - weight) * result.coefficients[i];
        }
      }
      return result;
    }

    static RealT distance(const HistoryPolynomial& a, const HistoryPolynomial& b)
    {
      RealT error = 0;
      for (size_t i = 0; i < N; ++i)
      {
        error = std::max(error, std::abs(a.coefficients[i] - b.coefficients[i]));
      }
      return error;
    }
  };

} // namespace GridKit::PhasorDynamics::detail
