#pragma once

#include <GridKit/Utilities/CliOptions/CliOptions.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    class CliOptionsTests
    {
    public:
      CliOptionsTests()
      {
      }

      virtual ~CliOptionsTests()
      {
      }

      /**
       * @brief Test for keys after parsing
       *
       * This method tests for specific keys after supplying a command-line
       * string to the CliOptions parser.
       */
      TestOutcome simpleParse()
      {
        using GridKit::Utilities::CliOptions;

        CommandLine cl{"app", "--file", "fakefile.txt", "-N", "16"};

        CliOptions options(cl.argc, cl.argv);

        TestStatus status{true};

        status *= options.getAppName() == "app";
        status *= options.hasKey("--file");
        status *= options.get<>("--file") == "fakefile.txt";
        status *= options.hasKey("-N");
        status *= options.get<int>("-N") == 16;

        status *= !options.hasKey("--bad");

        status *= throws<std::runtime_error>(
            [&]()
            { options.get<>("--bad"); });

        return status.report(__func__);
      }

    }; // class CliOptionsTests

  } // namespace Testing
} // namespace GridKit
