/**
 * @file Logger.hpp
 */

#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Manages and logs outputs from GridKit code.
     *
     * All methods and data in this class are static.
     *
     */
    class Logger
    {
    public:
      /// Enum specifying verbosity level for the output.
      enum Verbosity
      {
        NONE = 0,
        ERRORS,
        WARNINGS,
        SUMMARY,
        EVERYTHING
      };

      /// Buffer diagnostics on this thread and emit them as one record.
      /// Configure global output and verbosity before starting worker threads.
      class ScopedOutput
      {
      public:
        ScopedOutput();
        ~ScopedOutput();
        ScopedOutput(const ScopedOutput&)            = delete;
        ScopedOutput& operator=(const ScopedOutput&) = delete;

      private:
        std::ostringstream buffer_;
        std::ostream*      previous_;
      };

      // All methods and data are static so delete constructor and destructor.
      Logger()  = delete;
      ~Logger() = delete;

      static std::ostream& error();
      static std::ostream& warning();
      static std::ostream& summary();
      static std::ostream& misc();

      static void      setOutput(std::ostream& out);
      static void      openOutputFile(std::string filename);
      static void      closeOutputFile();
      static void      setVerbosity(Verbosity v);
      static void      raiseVerbosity(Verbosity v);
      static Verbosity verbosity();

      static std::vector<std::ostream*>& init();

    private:
      static std::ostream&              stream(Verbosity level);
      static thread_local std::ostream* thread_output_;
      static void                       updateVerbosity(std::vector<std::ostream*>& output_streams);

    private:
      static std::ostream               nullstream_;
      static std::ofstream              file_;
      static std::ostream*              logger_;
      static std::vector<std::ostream*> output_streams_;
      static std::vector<std::ostream*> tmp_;
      static Verbosity                  verbosity_;
    };
  } // namespace Utilities
} // namespace GridKit
