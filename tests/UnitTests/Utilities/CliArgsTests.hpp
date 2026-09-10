/**
 * @file CliArgsTests.hpp
 */

#pragma once

#include <sstream>
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

        CliArgs args{
          {.name = {"--opt1"}, .help = "an option with an argument"},

          {.name     = {"--opt2", "-o"},
           .help     = "An option with an argument and a short-hand",
           .type     = ArgType::Integer,
           .defaults = 0},

          {.name  = {"--opt3"},
           .help  = "An option with two arguments",
           .nargs = 2},

          {.name = {"--flag1", "-f"},
           .help = "A flag option (no argument)",
           .flag = true}
        };

        status *= args["opt1"].empty();
        status *= !args["opt2"].empty();
        status *= args["opt2"].as<int>() == 0;
        status *= args["opt3"].empty();
        status *= (&args["o"]) == (&args["opt2"]);
        status *= args["flag1"].as<bool>() == false;

        // bad: duplicate name
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Expect error because options cannot be duplicated\n";
        Log::setVerbosity(previous_verbosity);
        status *= throws<std::runtime_error>(
          [&]()
          { CliArgs args{{.name = {"--opt1"}}, {.name = {"--opt1", "-o"}}}; }
        );
        status *= throws<std::runtime_error>(
          [&]()
          {
            CliArgs args{
              {.name = {"--opt1", "-o"}}, {.name = {"--opt2", "-o"}}
            };
          }
        );

        return status.report(__func__);
      }

      /**
       * @brief Test for values after parsing
       */
      TestOutcome simpleParse()
      {
        using namespace GridKit::Utilities;

        CliArgs args{
          {.name = {"--file"}, .required = true},

          {.name = {"--twoargs"}, .nargs = 2},

          {.name = {"--size", "-N"}, .type = ArgType::Integer},

          {.name = {"--tol", "-t"}, .type = ArgType::Real},

          {.name     = {"--params", "-p"},
           .type     = ArgType::Real,
           .nargs    = 3,
           .defaults = {0, 0, 0}}
        };

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
          "-1.14",
          "-21"
        };

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
        status *= isEqual(y, -1.14);
        status *= isEqual(z, -21.0);
        status *=
          args.get<double, 3>("params") == args["params"].as<double, 3>();

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Expect error while testing that unrecognized options "
                       "are rejected\n";
        Log::setVerbosity(previous_verbosity);
        status *= throws<std::runtime_error>([&]()
                                             { args.get("bad"); });

        return status.report(__func__);
      }

      TestOutcome positionalAndRepeatable()
      {
        using namespace GridKit::Utilities;
        CliArgs args{
          {.name = {"study"}, .required = true, .positional = true},
          {.name = {"--item", "-i"}, .repeatable = true},
          {.name       = {"--pair"},
           .nargs      = 2,
           .defaults   = {0, 0},
           .repeatable = true},
          {.name = {"--flag", "-f"}, .flag = true}
        };
        CommandLine cl{
          "app",
          "-f",
          "--item",
          "friendly name=1",
          "study with spaces.json",
          "--pair",
          "-1",
          "2",
          "-i",
          "other=2",
          "--pair",
          "3",
          "4"
        };

        args.parseArgs(cl.argc, cl.argv);

        TestStatus status{true};
        status *= args.get("study") == "study with spaces.json";
        status *= args.get<bool>("flag");
        status *= args["item"].size() == 2;
        status *= args["item"][0]() == "friendly name=1";
        status *= args["item"][1]() == "other=2";
        status *= args["pair"].as<int, 4>() == std::array{-1, 2, 3, 4};
        std::ostringstream help;
        args.printHelp(help);
        status *= help.str().find("STUDY") != std::string::npos;
        status *= help.str().find("--study") == std::string::npos;
        status *= help.str().find("[-i ITEM]...") != std::string::npos;

        CommandLine escaped{"app", "--", "-study.json"};
        args.parseArgs(escaped.argc, escaped.argv);
        status *= args.get("study") == "-study.json";
        status *= !args.get<bool>("flag");
        status *= args["item"].empty();
        status *= args["pair"].as<int, 2>() == std::array{0, 0};

        CliArgs positionals{
          {.name = {"first"}, .positional = true},
          {.name = {"second"}, .nargs = 2, .positional = true},
          {.name = {"--flag"}, .flag = true}
        };
        CommandLine multiple{"app", "one", "--flag", "two", "three"};
        positionals.parseArgs(multiple.argc, multiple.argv);
        status *= positionals.get("first") == "one";
        status *= positionals["second"].as<2>() == std::array{"two"s, "three"s};
        status *= positionals.get<bool>("flag");
        return status.report(__func__);
      }

      TestOutcome positionalAndRepeatableErrors()
      {
        using namespace GridKit::Utilities;
        const auto argsAreBad = [&](auto&& cl)
        {
          CliArgs args{
            {.name = {"study"}, .required = true, .positional = true},
            {.name = {"--one"}},
            {.name = {"--item"}, .nargs = 2, .repeatable = true},
            {.name = {"--flag"}, .flag = true}
          };
          return throws<std::runtime_error>(
            [&]()
            { args.parseArgs(cl.argc, cl.argv); }
          );
        };
        TestStatus status{true};
        status *= argsAreBad(CommandLine{"app"});
        status *= argsAreBad(CommandLine{"app", "a.json", "b.json"});
        status *= argsAreBad(CommandLine{"app", "--study", "a.json"});
        status *=
          argsAreBad(CommandLine{"app", "a.json", "--one", "1", "--one", "2"});
        status *= argsAreBad(
          CommandLine{"app", "a.json", "--item", "1", "--item", "2", "3"}
        );
        status *= argsAreBad(
          CommandLine{"app", "a.json", "--item", "1", "2", "--item", "3"}
        );
        status *= argsAreBad(CommandLine{"app", "a.json", "--item", "--flag"});
        status *= argsAreBad(CommandLine{"app", "a.json", "--flag", "--flag"});
        status *= throws<std::runtime_error>(
          []
          {
            CliArgs args{{.name = {"pair"}, .nargs = 2, .positional = true}};

            CommandLine cl{"app", "one"};
            args.parseArgs(cl.argc, cl.argv);
          }
        );
        status *= throws<std::invalid_argument>(
          []
          {
            CliArgs args{
              {.name = {"file"}, .positional = true, .repeatable = true}
            };
          }
        );

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
          CliArgs args{
            {.name = {"--opt1"}, .required = true},

            {.name = {"--opt2", "-o"}, .nargs = 2},

            {.name = {"--flag1", "-f"}, .flag = true}
          };

          return throws<std::runtime_error>(
            [&]()
            { args.parseArgs(cl.argc, cl.argv); }
          );
        };

        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::EVERYTHING);

        // not providing required option
        Log::misc() << "Expect error for not providing required option\n";
        status *= argsAreBad(CommandLine{"app"});

        // giving an unexpected option
        Log::misc() << "Expect error about unexpected option\n";
        status *= argsAreBad(CommandLine{"app", "--opt1", "val", "-N", "16"});

        // incorrect number of values for an option
        Log::misc() << "Expect error because 1 < 2\n";
        status *=
          argsAreBad(CommandLine{"app", "--opt1", "val", "--opt2", "6"});
        Log::misc() << "Expect error because 3 > 2\n";
        status *= argsAreBad(
          CommandLine{"app", "--opt1", "val", "--opt2", "6", "7", "8"}
        );

        // passing a value to a flag option
        Log::misc() << "Expect error because flag options don't take args\n";
        status *= argsAreBad(CommandLine{"app", "--opt1", "val", "-f", "bad"});

        Log::misc() << "Expect error because option name can't be a number\n";
        status *= throws<std::runtime_error>(
          []
          { CliArgs args{{.name = {"--one", "-1"}}}; }
        );

        Log::setVerbosity(previous_verbosity);
        return status.report(__func__);
      }

    }; // class CliArgsTests
  } // namespace Testing
} // namespace GridKit
