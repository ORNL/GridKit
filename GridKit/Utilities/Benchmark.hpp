#pragma once

#include <chrono>
#include <format>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace GridKit
{
  namespace Utility
  {

    class Benchmark
    {
#ifdef GRIDKIT_BENCHMARK
      using Duration = std::chrono::duration<double, std::micro>;

      struct Run
      {
        struct Iteration
        {
          std::vector<Duration> observations_;
        };

        std::unordered_map<const char*, std::vector<Iteration>> iterations_;
        std::string                                             name_;
      };

      mutable std::mutex lock_;
      std::vector<Run>   runs_;

      Run& lastRun()
      {
        if (runs_.empty())
        {
          runs_.emplace_back();
        }

        return runs_.back();
      }

      Run::Iteration& lastIteration(const char* name)
      {
        auto& last_run = lastRun();

        auto iteration_iter = last_run.iterations_.find(name);
        if (iteration_iter == last_run.iterations_.end())
        {
          iteration_iter = std::get<0>(last_run.iterations_.insert({name, std::vector<Run::Iteration>(1)}));
        }

        std::vector<Run::Iteration>& iter_vec = std::get<1>(*iteration_iter);

        if (iter_vec.empty())
        {
          iter_vec.emplace_back();
        }

        return iter_vec.back();
      }

    public:
      void addTime(Duration duration, const char* name)
      {
        lock_.lock();

        auto& last_iter = lastIteration(name);
        last_iter.observations_.emplace_back(std::move(duration));

        lock_.unlock();
      }

      void newIteration()
      {
        lock_.lock();

        auto& last_run = lastRun();

        for (auto& var_pair : last_run.iterations_)
        {
          std::get<1>(var_pair).emplace_back();
        }

        lock_.unlock();
      }

      void newRun(std::string name)
      {
        if (runs_.size() > 0 && runs_.back().name_.empty())
        {
          runs_.back().name_ = std::move(name);
        }
        else
        {
          runs_.push_back({
              .iterations_ = {},
              .name_       = std::move(name),
          });
        }
      }

      void clear()
      {
        lock_.lock();

        runs_.clear();

        lock_.unlock();
      }

      std::string report() const
      {
        std::stringstream out;

        lock_.lock();

        for (const auto& run : runs_)
        {
          out << std::format("Run: {}\n", run.name_);
          constexpr std::string_view format_str = "{:>40} {:>15.4f} {:>15.4f} {:>6}\n";
          out << std::format("{:>40} {:>15} {:>15} {:>6}\n", "Variable", "Mean", "Std. Dev.", "N");

          for (const auto& var_pair : run.iterations_)
          {
            const auto& [var_name, iterations] = var_pair;

            double              sum                = 0;
            unsigned            total_observations = 0;
            std::vector<double> iter_sums;
            iter_sums.reserve(iterations.size());

            for (const auto& iter : iterations)
            {
              double iter_sum = 0;

              for (const auto& obs : iter.observations_)
              {
                sum      += obs.count();
                iter_sum += obs.count();

                total_observations++;
              }

              iter_sums.push_back(iter_sum);
            }

            double overall_mean = sum / total_observations;
            double iter_mean    = sum / iterations.size();
            double overall_var  = 0;
            double iter_var     = 0;

            for (size_t i = 0; i < iterations.size(); i++)
            {
              iter_var += std::pow(iter_sums[i] - iter_mean, 2) / (iterations.size() - 1);

              for (const auto& obs : iterations[i].observations_)
              {
                overall_var += std::pow(obs.count() - overall_mean, 2) / (total_observations - 1);
              }
            }

            out << std::format(format_str, var_name, overall_mean, std::sqrt(overall_var), total_observations);
            out << std::format(format_str, "-Total-", iter_mean, std::sqrt(iter_var), iterations.size());
          }
        }

        lock_.unlock();

        return out.str();
      }
#endif
    };

#ifdef GRIDKIT_BENCHMARK
    inline Benchmark benchmark;
#endif

    class Timer
    {
#ifdef GRIDKIT_BENCHMARK
      struct Inner
      {
        std::chrono::time_point<std::chrono::high_resolution_clock> start_;
        const char*                                                 name_;
      };

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
#ifdef GRIDKIT_BENCHMARK
        return Timer(std::chrono::high_resolution_clock::now(), name);
#else
        return Timer();
#endif
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
      benchmark.addTime(std::chrono::high_resolution_clock::now() - inner_timer->start_, inner_timer->name_);
#endif
    }
  }; // namespace Utility
}; // namespace GridKit