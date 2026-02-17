#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <vector>

#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    /**
     * @brief Represents a set of variables at a given time
     */
    struct OutputAtTime
    {
      double              t;
      std::vector<double> data;

      OutputAtTime(const std::vector<double>& v)
        : t{v[0]}
      {
        std::copy(next(begin(v)), end(v), std::back_inserter(data));
      }

      template <typename... TArgs>
      OutputAtTime(double t, TArgs&&... args)
        : t{t}, data{args...}
      {
      }

    private:
      friend OutputAtTime operator-(const OutputAtTime& a, const OutputAtTime& b)
      {
        assert(GridKit::Testing::isEqual(a.t, b.t, 1.0e-6));
        assert(a.data.size() == b.data.size());

        OutputAtTime ret(a);
        for (std::size_t i = 0; i < a.data.size(); ++i)
        {
          ret.data[i] -= b.data[i];
        }
        return ret;
      }

      friend double l1Norm(const OutputAtTime& a)
      {
        double ret = 0.0;
        for (auto v : a.data)
        {
          ret += std::abs(v);
        }
        return ret;
      }

      friend double l2Norm(const OutputAtTime& a)
      {
        double ret = 0.0;
        for (auto v : a.data)
        {
          ret += v * v;
        }
        return std::sqrt(ret);
      }

      friend double lInfNorm(const OutputAtTime& a)
      {
        double mx = 0.0;
        for (auto v : a.data)
        {
          mx = std::max(mx, std::abs(v));
        }
        return mx;
      }
    };

  } // namespace Testing
} // namespace GridKit
