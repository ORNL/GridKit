/**
 * @file CliArgsTests.hpp
 */

#pragma once

#include <string>

#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;
    using namespace std::string_literals;

    class CliArgsTests
    {
    public:
      CliArgsTests()
      {
      }

      virtual ~CliArgsTests()
      {
      }

      /**
       * @brief Test (good and bad) construction before parsing
       */
      TestOutcome construction()
      {
        using namespace GridKit::Utilities;

        TestStatus status{true};

        CliArgs args{{.name = {"--opt1"},
                      .help = "an option with an argument"},

                     {.name     = {"--opt2", "-o"},
                      .help     = "An option with an argument and a short-hand",
                      .type     = ArgType::Integer,
                      .defaults = 0},

                     {.name  = {"--opt3"},
                      .help  = "An option with two arguments",
                      .nargs = 2},

                     {.name = {"--flag1", "-f"},
                      .help = "A flag option (no argument)",
                      .flag = true}};

        status *= args["opt1"].empty();
        status *= !args["opt2"].empty();
        status *= args["opt2"].as<int>() == 0;
        status *= args["opt3"].empty();
        status *= (&args["o"]) == (&args["opt2"]);
        status *= args["flag1"].as<bool>() == false;

        // bad: duplicate name
        status *= throws<std::runtime_error>(
            [&]()
            {
              CliArgs args{
                  {.name = {"--opt1"}},
                  {.name = {"--opt1", "-o"}}};
            });
        status *= throws<std::runtime_error>(
            [&]()
            {
              CliArgs args{
                  {.name = {"--opt1", "-o"}},
                  {.name = {"--opt2", "-o"}}};
            });

        return status.report(__func__);
      }

      /**
       * @brief Test for values after parsing
       */
      TestOutcome simpleParse()
      {
        using namespace GridKit::Utilities;

        CliArgs args{{.name     = {"--file"},
                      .required = true},

                     {.name  = {"--twoargs"},
                      .nargs = 2},

                     {.name = {"--size", "-N"},
                      .type = ArgType::Integer},

                     {.name = {"--tol", "-t"},
                      .type = ArgType::Real},

                     {.name  = {"--params", "-p"},
                      .type  = ArgType::Real,
                      .nargs = 3}};

        CommandLine cl{
            "phony",
            "--file",
            "fake.txt",
            "--twoargs",
            "arg1",
            "arg2",
            "-N",
            "16",
            "-t",
            "0.01",
            "--params",
            "0.025",
            "1.14",
            "21"};

        args.parseArgs(cl.argc, cl.argv);

        TestStatus status{true};

        status *= args.getAppName() == "phony";
        status *= args["file"]() == "fake.txt";
        status *= args["twoargs"].as<2>() == std::array{"arg1"s, "arg2"s};

        auto [arg1, arg2] = args["twoargs"].as<2>();

        status *= arg1 == "arg1"s;
        status *= arg2 == "arg2"s;
        status *= args["N"].as<int>() == 16;

        auto [N] = args["N"].as<int, 1>();

        status *= N == 16;
        status *= args.get<int>("N") == 16;
        status *= args.get<double>("tol") == 0.01;

        auto [x, y, z] = args["params"].as<double, 3>();

        status *= isEqual(x, 0.025);
        status *= isEqual(y, 1.14);
        status *= isEqual(z, 21.0);
        status *=
            args.get<double, 3>("params") == args["params"].as<double, 3>();

        status *= throws<std::runtime_error>(
            [&]()
            { args.get("bad"); });

        return status.report(__func__);
      }

      /**
       * @brief Test for errors with bad command-line arguments
       */
      TestOutcome parsingErrors()
      {
        using namespace GridKit::Utilities;

        TestStatus status{true};

        auto argsAreBad = [&](auto&& cl)
        {
          CliArgs args{{.name     = {"--opt1"},
                        .required = true},

                       {.name  = {"--opt2", "-o"},
                        .nargs = 2},

                       {.name = {"--flag1", "-f"},
                        .flag = true}};

          return throws<std::runtime_error>(
              [&]()
              { args.parseArgs(cl.argc, cl.argv); });
        };

        Log::setVerbosity(Log::Verbosity::EVERYTHING);

        // not providing required option
        Log::misc() << "Expect error for not providing required option\n";
        status *= argsAreBad(CommandLine{"app"});

        // giving an unexpected option
        Log::misc() << "Expect error about unexpected option\n";
        status *= argsAreBad(CommandLine{"app", "--opt1", "val", "-N", "16"});

        // incorrect number of values for an option
        Log::misc() << "Expect error because 1 < 2\n";
        status *= argsAreBad(
            CommandLine{"app", "--opt1", "val", "--opt2", "6"});
        Log::misc() << "Expect error because 3 > 2\n";
        status *= argsAreBad(
            CommandLine{"app", "--opt1", "val", "--opt2", "6", "7", "8"});

        // passing a value to a flag option
        Log::misc() << "Expect error because flag options don't take args\n";
        status *= argsAreBad(CommandLine{"app", "--opt1", "val", "-f", "bad"});

        return status.report(__func__);
      }

    }; // class CliArgsTests

  } // namespace Testing
} // namespace GridKit
