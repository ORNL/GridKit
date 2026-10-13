#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
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
    };

    struct AffineSample
    {
      RealT constant;
      RealT weight;
    };

    struct Segment
    {
      RealT                begin;
      RealT                end;
      std::array<RealT, 4> coefficients;
      RealT                error{0};
      bool                 merge_failed{false}; ///< Unchanged pair already exceeds the merge budget.

      RealT evaluate(RealT t) const
      {
        auto        c = coefficients;
        const RealT x = (t - begin) / (end - begin);
        for (size_t n = 3; n > 0; --n)
        {
          for (size_t i = 0; i < n; ++i)
          {
            c[i] = (1 - x) * c[i] + x * c[i + 1];
          }
        }
        return c[0];
      }

      std::pair<Segment, Segment> split(RealT t) const
      {
        auto        c = coefficients;
        Segment     left{begin, t, {}}, right{t, end, {}};
        const RealT x         = (t - begin) / (end - begin);
        left.coefficients[0]  = c[0];
        right.coefficients[3] = c[3];
        for (size_t n = 3; n > 0; --n)
        {
          for (size_t i = 0; i < n; ++i)
          {
            c[i] = (1 - x) * c[i] + x * c[i + 1];
          }
          left.coefficients[4 - n]  = c[0];
          right.coefficients[n - 1] = c[n - 1];
        }
        return {left, right};
      }

      Segment restrictTo(RealT a, RealT b) const
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

      static RealT distance(const Segment& a, const Segment& b)
      {
        RealT error = 0;
        for (size_t i = 0; i < 4; ++i)
        {
          error = std::max(error, std::abs(a.coefficients[i] - b.coefficients[i]));
        }
        return error;
      }
    };

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
      else if (t != time())
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
      recent_.push_back({t, value});
    }

    /// Interpolation through accepted anchors and the current trial endpoint.
    static AffineSample trialSample(std::span<const Sample> past, RealT t, RealT q)
    {
      RealT a = 0;
      RealT b = 1;
      for (const auto& sample : past)
      {
        b *= (q - sample.time) / (t - sample.time);
      }
      for (size_t i = 0; i < past.size(); ++i)
      {
        RealT basis = (q - t) / (past[i].time - t);
        for (size_t j = 0; j < past.size(); ++j)
        {
          if (i != j)
          {
            basis *= (q - past[j].time) / (past[i].time - past[j].time);
          }
        }
        a += basis * past[i].value;
      }
      return {a, b};
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
      auto        value = [&](RealT q)
      {
        const auto a = trialSample(recent_, t, q);
        return a.constant + a.weight * endpoint;
      };
      const RealT c0 = recent_.back().value;
      const RealT c3 = endpoint;
      const RealT a  = 27 * value(begin + h / 3) - 8 * c0 - c3;
      const RealT b  = 27 * value(begin + 2 * h / 3) - c0 - 8 * c3;
      pending_       = Segment{begin, t, {c0, (2 * a - b) / 18, (2 * b - a) / 18, c3}};
      for (RealT coefficient : pending_->coefficients)
      {
        checkValue(coefficient);
      }
    }

    /// Compare both reconstruction and embedded interpolation at an interior point.
    RealT estimate(RealT q, RealT dense_value) const
    {
      checkValue(dense_value);
      const auto& p    = pending_.value();
      const RealT high = p.evaluate(q);
      checkValue(high);
      RealT error = std::abs(high - dense_value);
      if (recent_.size() > 1)
      {
        const auto low = trialSample(std::span<const Sample>(recent_).subspan(1), p.end, q);
        checkValue(low.constant + low.weight * p.coefficients[3]);
        error = std::max(error, std::abs(high - low.constant - low.weight * p.coefficients[3]));
      }
      // Reserve a compression budget only when merging is enabled.
      return error / errorScale(high, dense_value);
    }

    /// Check dense delayed output against the same history used by the residual.
    RealT estimateDelayed(RealT t, RealT delay, RealT output) const
    {
      checkValue(output);
      const auto& p        = pending_.value();
      const auto  a        = sample(p.end, p.end - (t - delay));
      const RealT expected = a.constant + a.weight * p.coefficients[3];
      return std::abs(output - expected) / errorScale(output, expected);
    }

    void accept()
    {
      const auto& p = pending_.value();
      segments_.push_back(p);
      recent_.push_back({p.end, p.coefficients[3]});
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
