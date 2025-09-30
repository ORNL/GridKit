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
        size_t                                                  current_iteration_;
      };

      mutable std::mutex lock_;
      std::vector<Run>   runs_;

      Run& lastRun()
      {
        if (runs_.empty())
        {
          runs_.emplace_back().current_iteration_ = 0;
        }

        return runs_.back();
      }

      Run::Iteration& lastIteration(const char* name)
      {
        auto& last_run = lastRun();

        auto iteration_iter = last_run.iterations_.find(name);
        if (iteration_iter == last_run.iterations_.end())
        {
          iteration_iter = std::get<0>(last_run.iterations_.insert({name, std::vector<Run::Iteration>()}));
        }

        std::vector<Run::Iteration>& iter_vec = std::get<1>(*iteration_iter);

        iter_vec.resize(last_run.current_iteration_);

        return iter_vec.back();
      }

    public:
      void addTime(Duration duration, const char* name)
      {
        const std::lock_guard guard(lock_);

        auto& last_iter = lastIteration(name);
        last_iter.observations_.emplace_back(std::move(duration));
      }

      void newIteration()
      {
        const std::lock_guard guard(lock_);

        auto& last_run = lastRun();

        for (auto& [var_name, iterations] : last_run.iterations_)
        {
          iterations.emplace_back();
        }

        last_run.current_iteration_++;
      }

      void newRun(std::string name)
      {
        const std::lock_guard guard(lock_);

        if (runs_.size() > 0 && runs_.back().name_.empty())
        {
          runs_.back().name_ = std::move(name);
        }
        else
        {
          runs_.push_back({
              .iterations_        = {},
              .name_              = std::move(name),
              .current_iteration_ = 0,
          });
        }
      }

      void clear()
      {
        const std::lock_guard guard(lock_);

        runs_.clear();
      }

      std::string report() const
      {
        std::stringstream out;

        const std::lock_guard guard(lock_);

        constexpr std::string_view format_str = "{:>40} {:>15.4f} {:>15.4f} {:>9}\n";
        constexpr std::string_view header_str = "{:-^40} {:-^15} {:-^15} {:-^9}\n";

        for (const auto& run : runs_)
        {
          out << std::format("Run: {}, ({} iterations)\n", run.name_, run.current_iteration_);
          out << std::format(header_str, "Variable", "Mean", "Std. Dev.", "N");

          std::unordered_map<std::string, std::vector<GridKit::Utility::Benchmark::Run::Iteration>> var_iterations;

          for (const auto& [var_name, iterations] : run.iterations_)
          {
            if (!var_iterations.contains(var_name))
            {
              var_iterations.insert({var_name, std::vector<GridKit::Utility::Benchmark::Run::Iteration>(run.current_iteration_)});
            }

            auto& new_iterations = var_iterations[var_name];

            for (size_t i = 0; i < iterations.size(); i++)
            {
              new_iterations[i].observations_.insert(new_iterations[i].observations_.end(), iterations[i].observations_.cbegin(), iterations[i].observations_.cend());
            }
          }

          for (const auto& [var_name, iterations] : var_iterations)
          {
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

            double overall_mean = sum / static_cast<double>(total_observations);
            double iter_mean    = sum / static_cast<double>(iterations.size());
            double overall_var  = 0;
            double iter_var     = 0;

            for (size_t i = 0; i < iterations.size(); i++)
            {
              iter_var += std::pow(iter_sums[i] - iter_mean, 2) / static_cast<double>(iterations.size() - 1);

              for (const auto& obs : iterations[i].observations_)
              {
                overall_var += std::pow(obs.count() - overall_mean, 2) / (total_observations - 1);
              }
            }

            out << std::format(format_str, var_name, overall_mean, std::sqrt(overall_var), total_observations);
            out << std::format(format_str, "-Total-", iter_mean, std::sqrt(iter_var), total_observations / iterations.size());
          }
        }

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

    template <bool ENABLE, typename F>
    auto time(const char* name, F func)
    {
      if constexpr (std::is_same<typename std::result_of<F&()>::type, void>::value)
      {
        auto timer = startTime<ENABLE>(name);
        func();
        endTime<ENABLE>(std::move(timer));
      }
      else
      {
        auto timer = startTime<ENABLE>(name);
        auto re    = func();
        endTime<ENABLE>(std::move(timer));

        return std::move(re);
      }
    }
  }; // namespace Utility
}; // namespace GridKit