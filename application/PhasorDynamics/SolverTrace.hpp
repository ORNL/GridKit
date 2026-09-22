#pragma once

#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <optional>
#include <vector>

#include <GridKit/Solver/Dynamic/Ida.hpp>

namespace GridKit::PhasorDynamics
{
  /// Optional segment-local IDA samples shared by the dynamics applications.
  class SolverTrace
  {
    using Stats = AnalysisManager::Sundials::IdaStats;

    struct Sample
    {
      std::size_t segment;
      const char* phase;
      double      t, h;
      Stats       stats;
    };

  public:
    explicit SolverTrace(const std::filesystem::path& path)
      : path_(path)
    {
    }

    template <class Solver>
    void record(const char* phase, double t, const Solver& solver, double h = 0.0)
    {
      if (!path_.empty())
      {
        samples_.push_back({segment_, phase, t, h, solver.getStats()});
      }
    }

    template <class Solver>
    void finish(double t, const Solver& solver)
    {
      record("end", t, solver);
      ++segment_;
    }

    template <class Solver>
    std::optional<std::function<void(double, double)>> callback(const Solver& solver)
    {
      if (path_.empty())
      {
        return {};
      }
      return [this, &solver](double t, double h)
      { record("step", t, solver, h); };
    }

    void write() const
    {
      if (path_.empty())
      {
        return;
      }
      if (!path_.parent_path().empty())
      {
        std::filesystem::create_directories(path_.parent_path());
      }
      std::ofstream output;
      output.exceptions(std::ios::failbit | std::ios::badbit);
      output.open(path_);
      output << "segment,phase,t,h,accepted_steps,residual_evals,jacobian_evals,error_test_failures\n"
             << std::setprecision(17);
      for (const auto& sample : samples_)
      {
        const auto& stats = sample.stats;
        output << sample.segment << ',' << sample.phase << ',' << sample.t << ',' << sample.h << ','
               << stats.num_steps_ << ',' << stats.num_residual_evals_ << ','
               << stats.num_jacobian_evals_ << ',' << stats.num_error_test_fails_ << '\n';
      }
      output.close();
    }

  private:
    std::filesystem::path path_;
    std::vector<Sample>   samples_;
    std::size_t           segment_ = 0;
  };
} // namespace GridKit::PhasorDynamics
