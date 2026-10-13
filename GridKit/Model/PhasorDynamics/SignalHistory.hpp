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

#include <GridKit/Model/PhasorDynamics/SignalHistoryPolynomial.hpp>

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

  public:
    using Cubic = detail::HistoryPolynomial<RealT, 4>;

    struct Segment : Cubic
    {
      RealT error{0};
      bool  merge_failed{false};
      bool  break_before{false}; ///< A known derivative corner cannot be merged away.

      Segment(RealT begin, RealT end, std::array<RealT, 4> coefficients)
        : Cubic{begin, end, coefficients}
      {
      }
    };

    /// Chebyshev-Lobatto nodes for degree-five polynomial reconstruction.
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

    void setPrehistory(std::function<RealT(RealT)> prehistory,
                       std::function<RealT(RealT)> derivative = {})
    {
      if (initialized())
      {
        throw std::logic_error("Set signal prehistory before initialization");
      }
      if (derivative && !prehistory)
      {
        throw std::invalid_argument("A prehistory derivative requires a value function");
      }
      prehistory_            = std::move(prehistory);
      prehistory_derivative_ = std::move(derivative);
    }

    /// Initial slope at the delayed query, not at the current source time.
    RealT initialDerivative(RealT t) const
    {
      if (!prehistory_)
      {
        return 0;
      }
      RealT derivative;
      if (prehistory_derivative_)
      {
        derivative = prehistory_derivative_(t);
      }
      else
      {
        // Backward differences stay inside the prescribed prehistory domain.
        const RealT dt = std::cbrt(std::numeric_limits<RealT>::epsilon())
                         * std::max(RealT{1}, std::abs(t));
        derivative = (3 * checkedPrehistory(t) - 4 * checkedPrehistory(t - dt)
                      + checkedPrehistory(t - 2 * dt))
                     / (2 * dt);
      }
      checkValue(derivative);
      return derivative;
    }

    /// Restart at an accepted continuous value after a known derivative corner.
    /// The first interval uses a secant; later intervals recover cubic order.
    void restartInterpolation()
    {
      if (!initialized() || pending_)
      {
        throw std::logic_error("Restart interpolation only at an accepted history time");
      }
      recent_[0]    = recent_[recent_count_ - 1];
      recent_count_ = 1;
      linear_start_ = true;
      trial_slope_.reset();
    }

    unsigned errorOrder() const
    {
      return linear_start_ ? 2 : static_cast<unsigned>(recent_count_ + 1);
    }

    bool hasPrehistory() const
    {
      return static_cast<bool>(prehistory_);
    }

    bool initialized() const
    {
      return recent_count_ != 0;
    }

    RealT time() const
    {
      return recent_[recent_count_ - 1].time;
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
      recent_count_       = 0;
      linear_start_       = false;
      compression_cursor_ = 0;
      trial_slope_.reset();
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
      initial_value_ = value;
      recent_[0]     = {t, value, initialDerivative(t)};
      recent_count_  = 1;
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
      return hermiteSample(past.back(), t, q, trialDerivative(past, t));
    }

    static AffineSample hermiteSample(const Sample& left, RealT t, RealT q, AffineSample derivative)
    {
      const RealT h   = t - left.time;
      const RealT x   = (q - left.time) / h;
      const RealT h00 = (1 + 2 * x) * (1 - x) * (1 - x);
      const RealT h10 = x * (1 - x) * (1 - x);
      const RealT h01 = x * x * (3 - 2 * x);
      const RealT h11 = x * x * (x - 1);
      return {h00 * left.value + h10 * h * left.derivative + h11 * h * derivative.constant,
              h01 + h11 * h * derivative.weight};
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
        return {recent_[recent_count_ - 1].value, 0};
      }
      if (q > time())
      {
        if (linear_start_)
        {
          const RealT weight = (q - time()) / (t - time());
          return {(1 - weight) * recent_[0].value, weight};
        }
        return hermiteSample(recent_[recent_count_ - 1], t, q, trialSlope(t));
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
      const RealT begin      = time();
      const RealT h          = t - begin;
      const auto  slope      = trialSlope(t);
      const RealT d1         = slope.constant + slope.weight * endpoint;
      const auto& left       = recent_[recent_count_ - 1];
      const RealT d0         = linear_start_ ? (endpoint - left.value) / h : left.derivative;
      pending_               = Segment{begin, t, {left.value, left.value + h * d0 / 3, endpoint - h * d1 / 3, endpoint}};
      pending_->break_before = linear_start_;
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

    /// Compare six polynomial samples against every crossed history segment.
    /// Prescribed source/prehistory functions still require adequate resolution.
    RealT estimateDelayed(RealT begin, RealT end, RealT delay, const std::array<RealT, 6>& values) const
    {
      for (RealT value : values)
      {
        checkValue(value);
      }
      const auto output = detail::HistoryPolynomial<RealT, 6>::interpolate(begin - delay, end - delay, values);
      return estimateDelayed(output);
    }

    /// Bound agreement with a polynomial over its overlaps with stored history.
    template <size_t N>
    RealT estimateDelayed(const detail::HistoryPolynomial<RealT, N>& output) const
    {
      static_assert(N >= 4);
      for (RealT coefficient : output.coefficients)
      {
        checkValue(coefficient);
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
        const auto cubic       = segment.restrictTo(a, b);
        const auto dense       = output.restrictTo(a, b);
        const auto expected    = cubic.template elevate<N>().coefficients;
        const auto [low, high] = std::minmax_element(expected.begin(), expected.end());
        const RealT minimum    = *low > 0 ? *low : (*high < 0 ? -*high : 0);
        for (size_t i = 0; i < N; ++i)
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
      if (recent_count_ == recent_.size())
      {
        std::move(recent_.begin() + 1, recent_.end(), recent_.begin());
        --recent_count_;
      }
      recent_[recent_count_++] = {p.end, p.coefficients[3], 3 * (p.coefficients[3] - p.coefficients[2]) / (p.end - p.begin)};
      linear_start_            = false;
      trial_slope_.reset();
      pending_.reset();
    }

    void discard() noexcept
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
        if (compression_cursor_ > 0)
        {
          --compression_cursor_;
        }
      }
      if (!options_.compress)
      {
        return;
      }

      auto protectedInterval = [&](RealT begin, RealT end)
      {
        if (end >= recent_[0].time)
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
      // Bound maintenance work per accepted step. Protected pairs are revisited
      // as their delayed queries move; failed pairs wait until a neighbor changes.
      const size_t work = std::min(size_t{32}, segments_.size());
      for (size_t attempt = 0; attempt < work && segments_.size() > 1; ++attempt)
      {
        if (compression_cursor_ + 1 >= segments_.size())
        {
          compression_cursor_ = 0;
        }
        size_t&     i     = compression_cursor_;
        const auto& left  = segments_[i];
        const auto& right = segments_[i + 1];
        if (left.merge_failed || right.break_before || protectedInterval(left.begin, right.end))
        {
          ++i;
          continue;
        }
        const RealT h = right.end - left.begin;
        Segment     merged{left.begin, right.end, {left.coefficients[0], left.coefficients[0] + h * (left.coefficients[1] - left.coefficients[0]) / (left.end - left.begin), right.coefficients[3] - h * (right.coefficients[3] - right.coefficients[2]) / (right.end - right.begin), right.coefficients[3]}};
        merged.break_before = left.break_before;
        merged.error        = std::max(left.error + Segment::distance(merged.restrictTo(left.begin, left.end), left),
                                right.error + Segment::distance(merged.restrictTo(right.begin, right.end), right));
        if (merged.error <= RealT{0.25} * options_.abs_tol)
        {
          segments_[i] = merged;
          segments_.erase(segments_.begin() + static_cast<std::ptrdiff_t>(i + 1));
          if (i > 0)
          {
            segments_[i - 1].merge_failed = false; // Its right neighbor changed.
            --i;
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
    AffineSample trialSlope(RealT t) const
    {
      if (!trial_slope_ || trial_slope_->first != t)
      {
        trial_slope_ = {t, trialDerivative(std::span(recent_.data(), recent_count_), t)};
      }
      return trial_slope_->second;
    }

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

    std::deque<Segment>                                   segments_;
    std::array<Sample, 3>                                 recent_{};
    size_t                                                recent_count_{0};
    size_t                                                compression_cursor_{0};
    bool                                                  linear_start_{false};
    std::optional<Segment>                                pending_;
    mutable std::optional<std::pair<RealT, AffineSample>> trial_slope_;
    std::function<RealT(RealT)>                           prehistory_;
    std::function<RealT(RealT)>                           prehistory_derivative_;
    RealT                                                 initial_time_{0};
    RealT                                                 initial_value_{0};
    SignalHistoryOptions<RealT>                           options_;
  };
} // namespace GridKit::PhasorDynamics
