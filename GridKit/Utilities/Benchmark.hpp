#pragma once

#include <chrono>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace GridKit
{
  namespace Utility
  {
    class Benchmark
    {
      using Duration = std::chrono::duration<double, std::micro>;

      struct Run
      {
        std::unordered_map<const char*, std::vector<Duration>> observations_;
      };

      std::mutex       lock_;
      std::vector<Run> runs_;

    public:
      void addTime(Duration duration, const char* name)
      {
        lock_.lock();

        if (runs_.empty())
        {
          runs_.emplace_back();
        }

        auto& last_run = runs_.back();

        if (auto found = last_run.observations_.find(name); found != last_run.observations_.end())
        {
          found->second.emplace_back(std::move(duration));
        }
        else
        {
          last_run.observations_.emplace(std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(1, std::move(duration)));
        }

        lock_.unlock();
      }

      void clear()
      {
        lock_.lock();

        runs_.clear();

        lock_.unlock();
      }
    };

    inline Benchmark benchmark;

    class Timer
    {
      struct Inner
      {
        std::chrono::time_point<std::chrono::high_resolution_clock> start_;
        const char*                                                 name_;
      };

#ifdef GRIDKIT_BENCHMARK
      std::unique_ptr<Inner> inner_;
#endif

    public:
      Timer()
      {
      }

#ifdef GRIDKIT_BENCHMARK
      Timer(std::chrono::time_point<std::chrono::high_resolution_clock>&& start, const char* name) : inner_(std::make_unique<Inner>(start, name))
      {
      }

      std::unique_ptr<Inner> intoInner() &&
      {
        return std::move(inner_);
      }
#else
      Timer(std::chrono::time_point<std::chrono::high_resolution_clock>&& start, const char* name)
      {
      }
#endif
    };

    template <bool ENABLE>
    Timer startTime(const char* name)
    {
      if constexpr (ENABLE)
      {
        return Timer(std::chrono::high_resolution_clock::now(), name);
      }
      else
      {
        return Timer();
      }
    }

    template <bool ENABLE>
    void endTime(Timer&& timer)
    {
      if constexpr (!ENABLE)
      {
        return;
      }

#ifdef GRIDKIT_BENCHMARK
      auto inner_timer = std::move(timer).intoInner();
      benchmark.addTime(std::chrono::high_resolution_clock::now() - inner_timer->start_, inner_timer->name);
#endif
    }
  }; // namespace Utility
}; // namespace GridKit