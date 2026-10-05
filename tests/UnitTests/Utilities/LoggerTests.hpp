/**
 * @file LoggerTests.hpp
 * @brief Contains definition of LoggerTests class.
 * @author Slaven Peles <peless@ornl.org>
 */

#pragma once
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <GridKit/Definitions.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Colors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Class implementing unit tests for Logger class.
     *
     * The LoggerTests class is implemented entirely in this header file.
     * Adding new unit test requires simply adding another method to this
     * class.
     */
    class LoggerTests
    {
    public:
      LoggerTests()
      {
      }

      virtual ~LoggerTests()
      {
      }

      /**
       * @brief Test the verbosity the Logger starts with.
       *
       * Developer mode (GRIDKIT_ENABLE_DEVELOPER_MODE) starts at EVERYTHING;
       * otherwise the Logger starts at WARNINGS. Must run before any test that
       * changes the verbosity.
       */
      TestOutcome defaultVerbosity()
      {
        using Log = GridKit::Utilities::Logger;

        TestStatus status;

#ifdef GRIDKIT_ENABLE_DEVELOPER_MODE
        status = (Log::verbosity() == Log::EVERYTHING);
#else
        status = (Log::verbosity() == Log::WARNINGS);
#endif

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for error log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams,
       * however only messages streamed to the error stream should be logged.
       */
      TestOutcome errorOutput()
      {
        using Log = GridKit::Utilities::Logger;
        std::string s1("Test error output ...");
        std::string s2("Another error output test ...\n");
        std::string answer = error_text() + s1 + "\n" + error_text() + s2;

        TestStatus status;

        std::ostringstream file;

        Log::setOutput(file);
        Log::setVerbosity(Log::ERRORS);
        Log::error() << s1 << std::endl;
        Log::error() << s2;

        Log::warning() << s1;
        Log::warning() << s2;
        Log::summary() << s1;
        Log::misc() << s1;

        status = (answer == file.str());

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for warning log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams,
       * however only messages streamed to the error and warning streams should
       * be logged.
       */
      TestOutcome warningOutput()
      {
        using Log = GridKit::Utilities::Logger;
        std::string s1("Test error output ...\n");
        std::string s2("Test warning output ...\n");
        std::string answer = error_text() + s1 + warning_text() + s2;

        TestStatus status;

        std::ostringstream file;

        Log::setOutput(file);
        Log::setVerbosity(Log::WARNINGS);

        Log::error() << s1;
        Log::warning() << s2;
        Log::summary() << s1;
        Log::misc() << s1;

        status = (answer == file.str());

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for result summary log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams,
       * however only messages streamed to the error, warning, and result summary
       * streams should be logged.
       */
      TestOutcome summaryOutput()
      {
        using Log = GridKit::Utilities::Logger;
        std::string s1("Test error output ...\n");
        std::string s2("Test warning output ...\n");
        std::string s3("Test summary output ...\n");
        std::string answer = error_text() + s1 + warning_text() + s2 + summary_ + s3;

        TestStatus status;

        std::ostringstream file;

        Log::setOutput(file);
        Log::setVerbosity(Log::SUMMARY);

        Log::error() << s1;
        Log::warning() << s2;
        Log::summary() << s3;
        Log::misc() << s1;

        status = (answer == file.str());

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for all other log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams
       * and all messages should be logged.
       */
      TestOutcome miscOutput()
      {
        using Log = GridKit::Utilities::Logger;
        std::string s1("Test error output ...\n");
        std::string s2("Test warning output ...\n");
        std::string s3("Test summary output ...\n");
        std::string s4("Test any other output ...\n");
        std::string answer = error_text() + s1 + warning_text() + s2 + summary_ + s3 + message_ + s4;

        TestStatus status;

        std::ostringstream file;

        Log::setOutput(file);
        Log::setVerbosity(Log::EVERYTHING);

        Log::error() << s1;
        Log::warning() << s2;
        Log::summary() << s3;
        Log::misc() << s4;

        status = (answer == file.str());

        return status.report(__func__);
      }

      /**
       * @brief Test raising the verbosity level.
       *
       * `raiseVerbosity` increases a lower verbosity to the requested level
       * and leaves an equal or higher verbosity unchanged.
       */
      TestOutcome raiseVerbosity()
      {
        using Log = GridKit::Utilities::Logger;

        TestStatus status;

        const auto previous_verbosity = Log::verbosity();

        Log::setVerbosity(Log::WARNINGS);
        Log::raiseVerbosity(Log::SUMMARY);
        status *= (Log::verbosity() == Log::SUMMARY);

        Log::raiseVerbosity(Log::ERRORS);
        status *= (Log::verbosity() == Log::SUMMARY);

        Log::setVerbosity(Log::EVERYTHING);
        Log::raiseVerbosity(Log::SUMMARY);
        status *= (Log::verbosity() == Log::EVERYTHING);

        Log::setVerbosity(previous_verbosity);

        return status.report(__func__);
      }

    private:
      /// Private method to return the string preceding error output
      std::string error_text()
      {
        using namespace Utilities::Colors;
        std::ostringstream stream;
        stream << "[" << RED << "ERROR" << CLEAR << "] ";
        return stream.str();
      }

      /// Private method to return the string preceding warning output
      std::string warning_text()
      {
        using namespace Utilities::Colors;
        std::ostringstream stream;
        stream << "[" << YELLOW << "WARNING" << CLEAR << "] ";
        return stream.str();
      }

      /// String preceding output of a result summary
      const std::string summary_ = "[SUMMARY] ";

      /// String preceding miscellaneous output
      const std::string message_ = "[MESSAGE] ";
    }; // class LoggerTests

  } // namespace Testing
} // namespace GridKit
