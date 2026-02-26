/**
 * @file CliArgs.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Define and parse expected command-line arguments
 */

#pragma once

#include <initializer_list>
#include <iostream>
#include <memory>
#include <string>

#include <GridKit/Utilities/CliArgs/ArgVector.hpp>
#include <GridKit/Utilities/CliArgs/Option.hpp>

namespace GridKit
{
  namespace Utilities
  {

    /**
     * @brief Table with which to specify and parse command-line arguments
     *
     * An application developer can specify expected options using a terse
     * syntax (using the keyword aggregate construction for each Option):
     *
     * ```
     * CliArgs args{{.name     = {"--opt1"},
     *               .help     = "an option with an argument"}
     *
     *              {.name     = {"--opt2", "-o"},
     *               .help     = "An option with an argument (and a short-hand)",
     *               .type     = ArgType::Int,
     *               .defaults = 0},
     *
     *              {.name     = {"--flag1", "-f"},
     *               .help     = "A flag option (no argument)",
     *               .flag     = true}};
     * ```
     *
     * @note The `{"--help", "-h"}` argument is added automatically. When it is
     * parsed from the command line. The help description will be printed and
     * the application will exit.
     *
     * @note The usage syntax is intended to be similar to
     * [ArgParse.jl](https://carlobaldassi.github.io/ArgParse.jl/stable/).
     *
     * @todo Allow option grouping for a more expressive help description.
     */
    class CliArgs
    {
    public:
      /**
       * @brief Construct with braced list (see class description)
       */
      CliArgs(std::initializer_list<Option> args);

      /// Destructor
      ~CliArgs();

      /**
       * @brief Parse command-line arguments using existing specification
       */
      void parseArgs(int argc, const char* argv[]);

      /**
       * @brief Print one-line application usage expression
       */
      void printUsage(std::ostream& os = std::cout) const;

      /**
       * @brief Print help description. Includes one-line usage expression and a
       * list of all expected options.
       */
      void printHelp(std::ostream& os = std::cout) const;

      /**
       * @brief Get the parsed values (or defaults) for the given argument.
       *
       * Arguments may be looked up with either their long-hand or their
       * short-hand forms, but without the hyphens (e.g., args["arg1"] using the
       * example from the class description).
       */
      const ArgVector& operator[](const std::string& name) const;

    private:
      /**
       * @brief Stream insertion operator for printing the argument table
       *
       * @note This is useful mainly for debugging purposes. It prints both the
       * list of options and the option name map, so the same Option will be
       * shown multiple times. Each Option is printed with its associated parsed
       * argument values.
       */
      friend std::ostream& operator<<(std::ostream& os, const CliArgs& args);

      std::unique_ptr<CliArgsImpl> pImpl_;
    };

  } // namespace Utilities
} // namespace GridKit
