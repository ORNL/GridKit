#pragma once

#include <algorithm>
#include <deque>
#include <functional>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename real_type>
    class SignalHistory
    {
    public:
      using RealT        = real_type;
      using PrehistoryFn = std::function<RealT(RealT)>;

      void requireWindow(RealT window)
      {
        if (window > window_)
        {
          window_ = window;
        }
      }

      RealT window() const
      {
        return window_;
      }

      bool active() const
      {
        return window_ > 0.0;
      }

      void reset()
      {
        samples_.clear();
      }

      void record(RealT time, RealT value)
      {
        if (!active())
        {
          return;
        }

        if (!samples_.empty() && time < samples_.back().time)
        {
          samples_.clear();
        }

        if (!samples_.empty() && time == samples_.back().time)
        {
          samples_.back().value = value;
          prune(time);
          return;
        }

        samples_.push_back(Sample{time, value});
        prune(time);
      }

      RealT read(RealT lookup_time, const PrehistoryFn& prehistory, RealT default_prehistory) const
      {
        if (samples_.empty())
        {
          return prehistory ? prehistory(lookup_time) : default_prehistory;
        }

        if (lookup_time < samples_.front().time)
        {
          return prehistory ? prehistory(lookup_time) : default_prehistory;
        }

        if (lookup_time == samples_.front().time)
        {
          return samples_.front().value;
        }

        if (lookup_time >= samples_.back().time)
        {
          return samples_.back().value;
        }

        auto upper = std::lower_bound(samples_.begin(),
                                      samples_.end(),
                                      lookup_time,
                                      [](const Sample& sample, RealT time)
                                      {
                                        return sample.time < time;
                                      });

        if (upper != samples_.end() && upper->time == lookup_time)
        {
          return upper->value;
        }

        auto       lower = std::prev(upper);
        const auto span  = upper->time - lower->time;
        if (span == 0.0)
        {
          return lower->value;
        }

        const auto theta = (lookup_time - lower->time) / span;
        return lower->value + theta * (upper->value - lower->value);
      }

    private:
      struct Sample
      {
        RealT time{0.0};
        RealT value{0.0};
      };

      void prune(RealT current_time)
      {
        const auto earliest_time = current_time - window_;
        while (samples_.size() > 2 && samples_[1].time <= earliest_time)
        {
          samples_.pop_front();
        }
      }

      std::deque<Sample> samples_;
      RealT              window_{0.0};
    };
  } // namespace PhasorDynamics
} // namespace GridKit
