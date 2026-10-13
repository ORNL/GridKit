#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace GridKit::PhasorDynamics
{
  template <class RealT>
  struct SignalHistoryOptions
  {
    RealT rel_tol{1e-5};
    RealT abs_tol{1e-7};
    bool  compress{false};
  };

  /// Accepted source history and affine interpolation within a trial step.
  /// Values in this class are real constants, never active AD variables.
  template <class RealT>
  class SignalHistory
  {
  public:
    struct Sample
    {
      RealT time;
      RealT value;
      RealT derivative{0};
    };

    struct AffineSample
    {
      RealT constant;
      RealT weight;
    };

  private:
    template <size_t N>
    struct Polynomial
    {
      RealT                begin;
      RealT                end;
      std::array<RealT, N> coefficients;
      RealT                error{0};
      bool                 merge_failed{false}; ///< Unchanged pair already exceeds the merge budget.

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

      std::pair<Polynomial, Polynomial> split(RealT t) const
      {
        auto        c = coefficients;
        Polynomial  left{begin, t, {}}, right{t, end, {}};
        const RealT x             = (t - begin) / (end - begin);
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

      Polynomial restrictTo(RealT a, RealT b) const
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

      static RealT distance(const Polynomial& a, const Polynomial& b)
      {
        RealT error = 0;
        for (size_t i = 0; i < N; ++i)
        {
          error = std::max(error, std::abs(a.coefficients[i] - b.coefficients[i]));
        }
        return error;
      }
    };

  public:
    using Segment = Polynomial<4>;

    /// Chebyshev-Lobatto nodes for reconstructing IDAS degree-five dense output.
    static constexpr std::array<RealT, 6> dense_nodes{0, 0.095491502812526288, 0.34549150281252627, 0.65450849718747373, 0.90450849718747373, 1};

    void setOptions(SignalHistoryOptions<RealT> options)
    {
      if (initialized())
      {
        throw std::logic_error("Set signal history options before initialization");
      }
      if (!std::isfinite(options.rel_tol) || options.rel_tol < 0 || !std::isfinite(options.abs_tol) || options.abs_tol <= 0)
      {
        throw std::invalid_argument("Signal history requires finite rtol >= 0 and atol > 0");
      }
      options_ = options;
    }

    void setPrehistory(std::function<RealT(RealT)> prehistory)
    {
      if (initialized())
      {
        throw std::logic_error("Set signal prehistory before initialization");
      }
      prehistory_ = std::move(prehistory);
    }

    bool hasPrehistory() const
    {
      return static_cast<bool>(prehistory_);
    }

    bool initialized() const
    {
      return !recent_.empty();
    }

    RealT time() const
    {
      return recent_.back().time;
    }

    size_t segmentCount() const
    {
      return segments_.size();
    }

    const std::deque<Segment>& segments() const
    {
      return segments_;
    }

    /// A fresh model initialization resets history but retains its specification.
    void reset()
    {
      recent_.clear();
      segments_.clear();
      pending_.reset();
    }

    void begin(RealT t)
    {
      if (!initialized())
      {
        initial_time_ = t;
      }
      // IDAS rounds a stop to its requested time within 100 unit roundoffs.
      else if (std::abs(t - time()) > 100 * std::numeric_limits<RealT>::epsilon()
                                          * std::max({RealT{1}, std::abs(t), std::abs(time())}))
      {
        throw std::logic_error("Signal history can only continue at its accepted time");
      }
    }

    void initialize(RealT t, RealT value)
    {
      begin(t);
      if (initialized())
      {
        return;
      }
      checkValue(value);
      if (prehistory_)
      {
        const RealT initial = checkedPrehistory(t);
        if (std::abs(initial - value) > options_.abs_tol + options_.rel_tol * std::max(std::abs(initial), std::abs(value)))
        {
          throw std::invalid_argument("Signal initial value is inconsistent with its prehistory");
        }
      }
      initial_value_   = value;
      RealT derivative = 0;
      if (prehistory_)
      {
        // A backward difference stays in the user-specified prehistory domain.
        const RealT dt = std::cbrt(std::numeric_limits<RealT>::epsilon()) * std::max(RealT{1}, std::abs(t));
        derivative     = (3 * value - 4 * checkedPrehistory(t - dt) + checkedPrehistory(t - 2 * dt)) / (2 * dt);
      }
      checkValue(derivative);
      recent_.push_back({t, value, derivative});
    }

    /// Endpoint derivative of the interpolant through recent values and U(t).
    static AffineSample trialDerivative(std::span<const Sample> past, RealT t)
    {
      RealT a = 0, b = 0;
      for (size_t i = 0; i < past.size(); ++i)
      {
        b           += 1 / (t - past[i].time);
        RealT basis  = 1 / (past[i].time - t);
        for (size_t j = 0; j < past.size(); ++j)
        {
          if (i != j)
          {
            basis *= (t - past[j].time) / (past[i].time - past[j].time);
          }
        }
        a += basis * past[i].value;
      }
      return {a, b};
    }

    /// C1 Hermite continuation, affine in the current source unknown U(t).
    static AffineSample trialSample(std::span<const Sample> past, RealT t, RealT q)
    {
      const auto& left = past.back();
      const RealT h    = t - left.time;
      const RealT x    = (q - left.time) / h;
      const RealT h00  = (1 + 2 * x) * (1 - x) * (1 - x);
      const RealT h10  = x * (1 - x) * (1 - x);
      const RealT h01  = x * x * (3 - 2 * x);
      const RealT h11  = x * x * (x - 1);
      const auto  d    = trialDerivative(past, t);
      return {h00 * left.value + h10 * h * left.derivative + h11 * h * d.constant,
              h01 + h11 * h * d.weight};
    }

    AffineSample sample(RealT t, RealT delay) const
    {
      const RealT q = t - delay;
      if (!initialized())
      {
        if (t != initial_time_)
        {
          throw std::logic_error("Configure signal delay step control before advancing the model");
        }
        if (prehistory_)
        {
          return {checkedPrehistory(q), 0};
        }
        return {0, 1}; // Equilibrium consistency solve: z = u.
      }
      if (q < initial_time_)
      {
        return {prehistory_ ? checkedPrehistory(q) : initial_value_, 0};
      }
      if (q == time())
      {
        return {recent_.back().value, 0};
      }
      if (q > time())
      {
        return trialSample(recent_, t, q);
      }
      if (q == initial_time_)
      {
        return {initial_value_, 0};
      }

      const auto it = std::lower_bound(segments_.begin(), segments_.end(), q, [](const Segment& s, RealT value)
                                       { return s.end < value; });
      if (it == segments_.end() || q < it->begin)
      {
        throw std::out_of_range("Delayed signal query is outside retained history");
      }
      return {it->evaluate(q), 0};
    }

    void prepare(RealT t, RealT endpoint)
    {
      if (!initialized() || !(t > time()))
      {
        throw std::logic_error("Signal history requires an advancing trial step");
      }
      checkValue(endpoint);
      const RealT begin = time();
      const RealT h     = t - begin;
      const auto  slope = trialDerivative(recent_, t);
      const RealT d1    = slope.constant + slope.weight * endpoint;
      const auto& left  = recent_.back();
      pending_          = Segment{begin, t, {left.value, left.value + h * left.derivative / 3, endpoint - h * d1 / 3, endpoint}};
      for (RealT coefficient : pending_->coefficients)
      {
        checkValue(coefficient);
      }
    }

    /// Compare the stored reconstruction with the candidate source trajectory.
    RealT estimate(RealT q, RealT dense_value) const
    {
      checkValue(dense_value);
      const auto& p    = pending_.value();
      const RealT high = p.evaluate(q);
      checkValue(high);
      const RealT error = std::abs(high - dense_value);
      // Reserve a compression budget only when merging is enabled.
      return error / errorScale(high, dense_value);
    }

    /// IDAS dense output has degree at most five. Bound its difference from
    /// every crossed cubic segment using Bernstein coefficients, including a
    /// query inside the current trial step. No delayed feature can be skipped
    /// merely because it lies between the six dense-output sample times.
    RealT estimateDelayed(RealT begin, RealT end, RealT delay, const std::array<RealT, 6>& values) const
    {
      // Inverse Bernstein collocation matrix: V_ij = binomial(5,j) x_i^j (1-x_i)^(5-j).
      static constexpr RealT inverse[6][6] = {
          {1, 0, 0, 0, 0, 0},
          {-2.4, 4.1888543819998318, -1.1577708763999663, 0.6111456180001682, -0.44222912360003364, 0.2},
          {2.52, -5.2944271909999161, 5.8499689437998486, -3.5055728090000842, 2.6300310562001514, -1.2},
          {-1.2, 2.6300310562001514, -3.5055728090000842, 5.8499689437998486, -5.2944271909999161, 2.52},
          {0.2, -0.44222912360003364, 0.6111456180001682, -1.1577708763999663, 4.1888543819998318, -2.4},
          {0, 0, 0, 0, 0, 1}};
      Polynomial<6> output{begin - delay, end - delay, {}};
      for (size_t i = 0; i < 6; ++i)
      {
        checkValue(values[i]);
        for (size_t j = 0; j < 6; ++j)
        {
          output.coefficients[i] += inverse[i][j] * values[j];
        }
        checkValue(output.coefficients[i]);
      }
      RealT error   = 0;
      auto  compare = [&](const Segment& segment)
      {
        const RealT a = std::max(output.begin, segment.begin);
        const RealT b = std::min(output.end, segment.end);
        if (!(b > a))
        {
          return;
        }
        const auto           cubic = segment.restrictTo(a, b);
        const auto           dense = output.restrictTo(a, b);
        // Elevate the cubic twice to compare polynomials in the same basis.
        std::array<RealT, 6> expected{};
        std::copy(cubic.coefficients.begin(), cubic.coefficients.end(), expected.begin());
        for (size_t degree = 4; degree <= 5; ++degree)
        {
          expected[degree] = expected[degree - 1];
          for (size_t i = degree - 1; i > 0; --i)
          {
            const RealT weight = static_cast<RealT>(i) / static_cast<RealT>(degree);
            expected[i]        = weight * expected[i - 1] + (1 - weight) * expected[i];
          }
        }
        const auto [low, high] = std::minmax_element(expected.begin(), expected.end());
        const RealT minimum    = *low > 0 ? *low : (*high < 0 ? -*high : 0);
        for (size_t i = 0; i < 6; ++i)
        {
          error = std::max(error, std::abs(dense.coefficients[i] - expected[i]) / errorScale(minimum, minimum));
        }
      };
      const auto first = std::upper_bound(segments_.begin(), segments_.end(), output.begin, [](RealT t, const Segment& segment)
                                          { return t < segment.end; });
      for (auto it = first; it != segments_.end() && it->begin < output.end; ++it)
      {
        compare(*it);
      }
      compare(pending_.value());
      if (output.begin < initial_time_)
      {
        // Arbitrary analytic prehistory is sampled, just like a source function.
        const RealT last = std::min(output.end, initial_time_);
        for (size_t i = 0; i <= 6; ++i)
        {
          const RealT q        = output.begin + (last - output.begin) * static_cast<RealT>(i) / 6;
          const RealT expected = prehistory_ ? checkedPrehistory(q) : initial_value_;
          const RealT actual   = output.evaluate(q);
          error                = std::max(error, std::abs(actual - expected) / errorScale(actual, expected));
        }
      }
      return error;
    }

    void accept()
    {
      const auto& p = pending_.value();
      segments_.push_back(p);
      recent_.push_back({p.end, p.coefficients[3], 3 * (p.coefficients[3] - p.coefficients[2]) / (p.end - p.begin)});
      if (recent_.size() > 3)
      {
        recent_.erase(recent_.begin());
      }
      pending_.reset();
    }

    void discard()
    {
      pending_.reset();
    }

    /// Called after dense output has been published for the accepted step.
    void maintain(std::span<const RealT> delays)
    {
      if (delays.empty() || !initialized())
      {
        return;
      }
      const RealT oldest = time() - *std::max_element(delays.begin(), delays.end());
      while (!segments_.empty() && segments_.front().end < oldest)
      {
        segments_.pop_front();
      }
      if (!options_.compress)
      {
        return;
      }

      auto protectedInterval = [&](RealT begin, RealT end)
      {
        if (end >= recent_.front().time)
        {
          return true;
        }
        for (RealT delay : delays)
        {
          if (begin <= time() - delay && time() - delay <= end)
          {
            return true;
          }
        }
        return false;
      };
      for (size_t i = 0; i + 1 < segments_.size();)
      {
        const auto& left  = segments_[i];
        const auto& right = segments_[i + 1];
        if (left.merge_failed || protectedInterval(left.begin, right.end))
        {
          ++i;
          continue;
        }
        const RealT h = right.end - left.begin;
        Segment     merged{left.begin, right.end, {left.coefficients[0], left.coefficients[0] + h * (left.coefficients[1] - left.coefficients[0]) / (left.end - left.begin), right.coefficients[3] - h * (right.coefficients[3] - right.coefficients[2]) / (right.end - right.begin), right.coefficients[3]}};
        merged.error = std::max(left.error + Segment::distance(merged.restrictTo(left.begin, left.end), left),
                                right.error + Segment::distance(merged.restrictTo(right.begin, right.end), right));
        if (merged.error <= RealT{0.25} * options_.abs_tol)
        {
          segments_[i] = merged;
          segments_.erase(segments_.begin() + static_cast<std::ptrdiff_t>(i + 1));
          if (i > 0)
          {
            segments_[i - 1].merge_failed = false; // Its right neighbor changed.
          }
        }
        else
        {
          segments_[i].merge_failed = true;
          ++i;
        }
      }
    }

  private:
    RealT errorScale(RealT a, RealT b) const
    {
      return (options_.compress ? RealT{0.75} : RealT{1}) * options_.abs_tol
             + options_.rel_tol * std::max(std::abs(a), std::abs(b));
    }

    static void checkValue(RealT value)
    {
      if (!std::isfinite(value))
      {
        throw std::runtime_error("Nonfinite delayed signal value");
      }
    }

    RealT checkedPrehistory(RealT t) const
    {
      const RealT value = prehistory_(t);
      checkValue(value);
      return value;
    }

    std::deque<Segment>         segments_;
    std::vector<Sample>         recent_;
    std::optional<Segment>      pending_;
    std::function<RealT(RealT)> prehistory_;
    RealT                       initial_time_{0};
    RealT                       initial_value_{0};
    SignalHistoryOptions<RealT> options_;
  };
} // namespace GridKit::PhasorDynamics
