/**
 * @file CliArgs.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Define and parse expected command-line arguments
 */

#pragma once

#include <array>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace GridKit
{
  namespace Utilities
  {
    struct CliArgsImpl;

    /**
     * @brief Represents a single value of arbitrary type
     *
     * The value is internally represented as a std::string that is parsed when
     * requested as a specific type
     */
    class ArgValue
    {
    public:
      /**
       * @brief Default construction results in empty value
       */
      ArgValue() = default;

      /**
       * @brief Construct from any value
       */
      template <typename T>
      ArgValue(T&& val)
        : value_((std::stringstream() << val).str())
      {
      }

      /**
       * @brief Assign any value
       */
      template <typename T>
      ArgValue& operator=(T&& val)
      {
        value_ = (std::stringstream() << val).str();
        return *this;
      }

      /**
       * @brief Check if no value is contained
       */
      bool empty() const
      {
        return value_.empty();
      }

      /**
       * @brief Get string representation
       */
      const std::string& get() const
      {
        return value_;
      }

      /**
       * @brief Get string representation
       */
      const std::string& operator()() const
      {
        return get();
      }

      /**
       * @brief Get value of specific expected type
       */
      template <typename T>
      T as() const
      {
        T ret;
        std::stringstream(value_) >> ret;
        return ret;
      }

    private:
      /// Internal representation
      std::string value_{};
    };

    namespace Detail
    {
      template <typename T>
      struct ArgValGetHelper
      {
        T operator()(const ArgValue& val) const
        {
          return val.as<T>();
        }
      };

      template <>
      struct ArgValGetHelper<std::string>
      {
        const std::string& operator()(const ArgValue& val) const
        {
          return val.get();
        }
      };
    } // namespace Detail

    /**
     * @brief Get internal value from an ArgValue object as a specified type.
     *
     * The implementation is specialized via Detail::ArgValGetHelper simply in
     * order to avoid a copy for the case when T is std::string
     */
    template <typename T>
    decltype(auto) get(const ArgValue& val)
    {
      return Detail::ArgValGetHelper<T>{}(val);
    }

    /**
     * @brief Represents a set of 0 or more values associated with a given
     * command-line option
     *
     * This class allows the internal data to be interpreted as a discrete set
     * of values (std::array) of a given type or as a single value, hopefully
     * making it more intuitive for users. For example...
     *
     * If an ArgVector v is expected to hold one value of type `double`, that
     * value can be extracted with `auto r = v.as<double>();`. This means `r`
     * will be a `double` value assigned with the contents of the first element
     * of `v`.
     *
     * If an ArgVector v is expected to hold two values of type `int`, on the
     * other hand, the values can be extracted as a `std::array<int, 2>` using
     * `as<int, 2>()`. Using `std::array` allows structured bindings for the
     * user: `auto [a, b] = v.as<int, 2>();`
     */
    class ArgVector
    {
    public:
      /**
       * @brief Default is empty (no values)
       */
      ArgVector() = default;

      /**
       * @brief Construct with a single value
       */
      template <typename T>
      ArgVector(T&& val)
        : vec{val}
      {
      }

      /**
       * @brief Interpret as `N` values of type `T`
       */
      template <typename T, std::size_t N>
      std::array<T, N> as() const
      {
        assert(vec.size() == N);
        std::array<T, N> ret;
        for (std::size_t i = 0; i < N; ++i)
        {
          ret[i] = vec[i].as<T>();
        }
        return ret;
      }

      /**
       * @brief Interpret as `N` values of type `std::string`
       */
      template <std::size_t N>
      decltype(auto) as() const
      {
        return as<std::string, N>();
      }

      /**
       * @brief Interpret as single value of type `T`
       */
      template <typename T>
      decltype(auto) as() const
      {
        return get<T>(vec[0]);
      }

      /**
       * @brief Interpret as single `std::string` value
       */
      const std::string& operator()() const
      {
        return get<std::string>(vec[0]);
      }

      /**
       * @brief Get one of the contained `ArgValue` objects
       */
      const ArgValue& operator[](std::size_t i) const
      {
        return vec[i];
      }

      /**
       * @brief Check for existence of values
       */
      bool empty() const
      {
        return vec.empty();
      }

      /**
       * @brief Get number of `ArgValue` objects
       */
      std::size_t size() const
      {
        return vec.size();
      }

    private:
      /// Internal set of values
      std::vector<ArgValue> vec;

      friend class CliArgsImpl;
    };

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
