/**
 * @file Option.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 */

#pragma once

#include <array>
#include <string>

#include <GridKit/Utilities/CliArgs/ArgVector.hpp>

namespace GridKit
{
  namespace Utilities
  {

    /**
     * @brief Specify expected type to which an argument should parse
     *
     * @note Currently this is only used for the "help" output to let users
     * know what type of value to provide.
     */
    enum class ArgType
    {
      String,     ///< Expect string value
      Real,       ///< Expect floating-point value
      Integer,    ///< Expect integer value
      Boolean,    ///< Expect bool value
      Unspecified ///< Since the `type` member is optional
    };

    /**
     * @brief Specifies a command-line option's parameters
     *
     * @note An open struct is used for this to enable keyword construction
     * (e.g.,
     * `Option opt{.name = {"--arg", "-a"}, .help = "Option description"};`
     * ). See CliArgs.
     *
     * @todo
     * - Allow positional arguments
     * - Allow `nargs` to be "arbitrary"
     */
    struct Option
    {
      /**
       * @brief Key(s) by which to refer to this Option
       *
       * This is the only field that must be explicitly initialized. It expects
       * 1 or 2 elements.  The first element must exist and must be the
       * long-hand form, `--name` (must have exactly two hyphens).  The second
       * element is optional but must be the short-hand form, `-n` (that is,
       * only one `-`).
       */
      std::array<std::string, 2> name;
      /// Help description to be printed on request. Can be multiple lines.
      std::string                help{};
      /// Arguments are optional be default. This makes an argument required
      bool                       required{false};
      /// The expected type to which a value should parse
      ArgType                    type{ArgType::Unspecified};
      /**
       * @brief This indicates that the option is a flag, and not followed by a
       * value
       *
       * Flag options are automatically configured to be false by default and
       * flip to true if parsed.
       */
      bool                       flag{false};
      /// The number of values that can be provided with an option
      std::size_t                nargs{1};
      /// Default value(s)
      ArgVector                  defaults{};
    };

  } // namespace Utilities
} // namespace GridKit
