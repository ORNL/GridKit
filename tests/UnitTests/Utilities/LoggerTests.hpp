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
        using out = GridKit::Utilities::Logger;

        TestStatus status;

#ifdef GRIDKIT_ENABLE_DEVELOPER_MODE
        status = (out::verbosity() == out::EVERYTHING);
#else
        status = (out::verbosity() == out::WARNINGS);
#endif

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for error log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams,
       * however only mesages streamed to the error stream should be logged.
       */
      TestOutcome errorOutput()
      {
        using out = GridKit::Utilities::Logger;
        std::string s1("Test error output ...");
        std::string s2("Another error output test ...\n");
        std::string answer = error_text() + s1 + "\n" + error_text() + s2;

        TestStatus status;

        std::ostringstream file;

        out::setOutput(file);
        out::setVerbosity(out::ERRORS);
        out::error() << s1 << std::endl;
        out::error() << s2;

        out::warning() << s1;
        out::warning() << s2;
        out::summary() << s1;
        out::misc() << s1;

        status = (answer == file.str());

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for warning log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams,
       * however only mesages streamed to the error and warning streams should
       * be logged.
       */
      TestOutcome warningOutput()
      {
        using out = GridKit::Utilities::Logger;
        std::string s1("Test error output ...\n");
        std::string s2("Test warning output ...\n");
        std::string answer = error_text() + s1 + warning_text() + s2;

        TestStatus status;

        std::ostringstream file;

        out::setOutput(file);
        out::setVerbosity(out::WARNINGS);

        out::error() << s1;
        out::warning() << s2;
        out::summary() << s1;
        out::misc() << s1;

        status = (answer == file.str());

        return status.report(__func__);
      }

      /**
       * @brief Test data stream for result summary log messages.
       *
       * This method tests streaming messages to `Logger::error()` data
       * stream. The method streams messages to all available output streams,
       * however only mesages streamed to the error, warning, and result summary
       * streams should be logged.
       */
      TestOutcome summaryOutput()
      {
        using out = GridKit::Utilities::Logger;
        std::string s1("Test error output ...\n");
        std::string s2("Test warning output ...\n");
        std::string s3("Test summary output ...\n");
        std::string answer = error_text() + s1 + warning_text() + s2 + summary_ + s3;

        TestStatus status;

        std::ostringstream file;

        out::setOutput(file);
        out::setVerbosity(out::SUMMARY);

        out::error() << s1;
        out::warning() << s2;
        out::summary() << s3;
        out::misc() << s1;

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
        using out = GridKit::Utilities::Logger;
        std::string s1("Test error output ...\n");
        std::string s2("Test warning output ...\n");
        std::string s3("Test summary output ...\n");
        std::string s4("Test any other output ...\n");
        std::string answer = error_text() + s1 + warning_text() + s2 + summary_ + s3 + message_ + s4;

        TestStatus status;

        std::ostringstream file;

        out::setOutput(file);
        out::setVerbosity(out::EVERYTHING);

        out::error() << s1;
        out::warning() << s2;
        out::summary() << s3;
        out::misc() << s4;

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
        using out = GridKit::Utilities::Logger;

        TestStatus status;

        const auto previous_verbosity = out::verbosity();

        out::setVerbosity(out::WARNINGS);
        out::raiseVerbosity(out::SUMMARY);
        status *= (out::verbosity() == out::SUMMARY);

        out::raiseVerbosity(out::ERRORS);
        status *= (out::verbosity() == out::SUMMARY);

        out::setVerbosity(out::EVERYTHING);
        out::raiseVerbosity(out::SUMMARY);
        status *= (out::verbosity() == out::EVERYTHING);

        out::setVerbosity(previous_verbosity);

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
