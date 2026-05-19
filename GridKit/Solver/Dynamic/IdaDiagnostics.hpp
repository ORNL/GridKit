#pragma once

#include <filesystem>
#include <optional>
#include <stdexcept>
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

    IdaStats idaStatsDelta(const IdaStats& end_stats, const IdaStats& start_stats);

    void writeIdaStatsJson(const IdaStatsReport& report, const IdaDiagnosticsOutput& output);
  } // namespace Sundials
} // namespace AnalysisManager
