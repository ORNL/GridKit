#pragma once

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <GridKit/Solver/Dynamic/Ida.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    struct IdaStatsSegment
    {
      double   start_time{};
      double   end_time{};
      int      output_steps{};
      IdaStats stats;
    };

    struct IdaDiagnosticsOutput
    {
      std::filesystem::path        file;
      std::optional<IdaLogOptions> log;
    };

    struct IdaStatsReport
    {
      IdaStats                     summary;
      std::vector<IdaStatsSegment> segments;
      std::optional<IdaLogOptions> log;
    };

    struct IdaStepSample
    {
      std::size_t segment_index{};
      std::size_t segment_step_index{};
      std::size_t global_step{};
      double      step_start_time{};
      double      step_end_time{};
      double      last_step{};
      double      current_step{};
      int         last_order{};
      int         current_order{};
      double      current_cj{};
      IdaStats    counter_delta;
    };

    struct IdaStepHistorySegment
    {
      std::size_t                segment_index{};
      double                     start_time{};
      double                     end_time{};
      int                        output_steps{};
      IdaStats                   final_stats;
      std::vector<IdaStepSample> steps;
    };

    struct IdaStepHistoryReport
    {
      std::vector<IdaStepHistorySegment> segments;
    };

    class IdaStatsRecorder
    {
    public:
      explicit IdaStatsRecorder(bool enabled);

      bool enabled() const;

      template <class IdaT>
      void beginSegment(const IdaT& ida)
      {
        if (!enabled_)
        {
          return;
        }
        segment_start_stats_ = ida.getStats();
      }

      template <class IdaT>
      void endSegment(const IdaT& ida, double start_time, double end_time, int output_steps)
      {
        if (!enabled_)
        {
          return;
        }
        if (!segment_start_stats_.has_value())
        {
          throw std::logic_error("IDA stats segment ended before it was started");
        }

        const auto end_stats = ida.getStats();
        recordSegment(*segment_start_stats_, end_stats, start_time, end_time, output_steps);
        segment_start_stats_.reset();
      }

      IdaStatsReport report(std::optional<IdaLogOptions> log = {}) const;

    private:
      void recordSegment(const IdaStats& start_stats,
                         const IdaStats& end_stats,
                         double          start_time,
                         double          end_time,
                         int             output_steps);

      bool                         enabled_{false};
      std::optional<IdaStats>      segment_start_stats_;
      IdaStats                     summary_;
      std::vector<IdaStatsSegment> segments_;
    };

    class IdaStepHistoryRecorder
    {
    public:
      explicit IdaStepHistoryRecorder(bool enabled);

      bool enabled() const;

      template <class IdaT>
      void beginSegment(const IdaT& ida, double start_time, double end_time, int output_steps)
      {
        if (!enabled_)
        {
          return;
        }

        current_segment_.emplace();
        current_segment_->segment_index = segments_.size();
        current_segment_->start_time    = start_time;
        current_segment_->end_time      = end_time;
        current_segment_->output_steps  = output_steps;
        previous_step_stats_            = ida.getStats();
      }

      void recordStep(const IdaStats& stats);

      template <class IdaT>
      void endSegment(const IdaT& ida)
      {
        if (!enabled_)
        {
          return;
        }
        if (!current_segment_.has_value() || !previous_step_stats_.has_value())
        {
          throw std::logic_error("IDA step-history segment ended before it was started");
        }

        current_segment_->final_stats = ida.getStats();
        segments_.push_back(std::move(*current_segment_));
        current_segment_.reset();
        previous_step_stats_.reset();
      }

      IdaStepHistoryReport report() const;

    private:
      bool                                 enabled_{false};
      std::optional<IdaStats>              previous_step_stats_;
      std::optional<IdaStepHistorySegment> current_segment_;
      std::vector<IdaStepHistorySegment>   segments_;
      std::size_t                          accepted_steps_{0};
    };

    IdaStats idaStatsDelta(const IdaStats& end_stats, const IdaStats& start_stats);

    void writeIdaStatsJson(const IdaStatsReport& report, const IdaDiagnosticsOutput& output);
    void writeIdaStepHistoryJson(const IdaStepHistoryReport& report, const IdaDiagnosticsOutput& output);
  } // namespace Sundials
} // namespace AnalysisManager
