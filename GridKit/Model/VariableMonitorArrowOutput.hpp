#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <string>

#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace Model
  {
    class VariableMonitorArrowOutputImpl;

    /**
     * @brief Monitor output writing Apache Arrow IPC
     *
     * Writes either the IPC file format (Feather v2) or the IPC stream
     * format. A file name is required; the name may refer to a FIFO to
     * stream record batches to a live consumer.
     *
     * The file format needs its footer, written in stop(), to be readable;
     * a run that dies before stopping the monitor leaves the file without
     * one. The stream format has no footer and stays readable up to the
     * last flushed batch.
     *
     * All Arrow-dependent code lives in the implementation, which is only
     * available when GridKit is built with `GridKit_ENABLE_ARROW=ON`;
     * otherwise construction throws.
     */
    class VariableMonitorArrowOutput
    {
    public:
      VariableMonitorArrowOutput() = delete;

      /**
       * @brief Construct output for the given output file
       *
       * @param file_name Output file name (may be a FIFO for live streaming)
       * @param stream_format Write the IPC stream format instead of the IPC
       * file format
       * @param layout Wide or long table layout
       * @param batch_rows Rows buffered before a record batch is written and
       * flushed
       */
      VariableMonitorArrowOutput(const std::string& file_name, bool stream_format, VariableMonitorLayout layout, int64_t batch_rows);

      VariableMonitorArrowOutput(const VariableMonitorArrowOutput&) = delete;

      VariableMonitorArrowOutput(VariableMonitorArrowOutput&&) noexcept;

      ~VariableMonitorArrowOutput();

      /// Open the output file in binary mode
      void start();

      /**
       * @brief Build the schema from the monitored variables and create the
       * IPC writer
       *
       * Const because it is called on const sinks while printing; only
       * implementation state changes (compare `Json::after_first`).
       */
      void beginTable(std::span<const VariableMonitorBase::Entry> entries) const;

      /// Buffer one sample, writing full record batches as needed
      void append(std::size_t sample, double t, std::span<const VariableMonitorBase::Entry> entries) const;

      /// Flush buffered rows, finalize the IPC output, and close the file
      void stop();

    private:
      // Members are only read by the Arrow-enabled implementation
      /// Output file name
      [[maybe_unused]] std::string file_name_;

      /// Write the IPC stream format instead of the IPC file format
      [[maybe_unused]] bool stream_format_{false};

      /// Table layout
      [[maybe_unused]] VariableMonitorLayout layout_{VariableMonitorLayout::WIDE};

      /// Rows buffered before a record batch is written and flushed
      [[maybe_unused]] int64_t batch_rows_{256};

      /// Arrow-dependent implementation; defined in the only TU that
      /// includes Arrow headers
      mutable std::unique_ptr<VariableMonitorArrowOutputImpl> impl_;
    };
  } // namespace Model
} // namespace GridKit
