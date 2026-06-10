#include <GridKit/Utilities/Stream.hpp>
#include <GridKit/Utilities/String.hpp>

namespace GridKit
{
  namespace Utilities
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Associate an Option with its parsed arguments
     */
    struct OptionArgumentPair
    {
      /**
       * @brief Construct with an Option spec and an empty set of values
       */
      explicit OptionArgumentPair(Option opt)
        : option(std::move(opt))
      {
      }

      // OptionArgumentPair(OptionArgumentPair

      /// Reference to existing option
      Option    option;
      /// Parsed values
      ArgVector values{};
    };

    /**
     * @brief Implementation of CliArgs
     */
    struct CliArgsImpl
    {
      ///@{
      /// @brief Implmentation of CliArgs
      CliArgsImpl(std::initializer_list<Option> args);

      void parseArgs(int argc, const char* argv[]);
      void printUsage(std::ostream& os = std::cout) const;
      void printHelp(std::ostream& os = std::cout) const;

      const ArgVector& operator[](const std::string& name) const;
      ///@}

      /// Map a name to an Option (fail if duplicate name)
      void mapName(const std::string& name, OptionArgumentPair* opt);
      /// Perform sanity assertions and add Option to name_map_
      void setupArg(OptionArgumentPair& opt);

      /// Option lookup via (hyphenless) name
      const OptionArgumentPair* findByName(const std::string& name) const;

      /// Option lookup via (hyphenless) name
      OptionArgumentPair* findByName(const std::string& name)
      {
        return const_cast<OptionArgumentPair*>(
            const_cast<const CliArgsImpl*>(this)->findByName(name));
      }

      /// Option lookup via raw (hyphened) name
      OptionArgumentPair* findByRawName(const std::string& raw);

      /// Store application name ("app" placeholder is changed in parseArgs())
      std::string                                          app_name_{"app"};
      /// Collection of options
      std::vector<OptionArgumentPair>                      table_;
      /// Map of names to options
      std::unordered_map<std::string, OptionArgumentPair*> name_map_;
    };

    /**
     * @brief Print all the values in an ArgVector to a string
     *
     * @note This is useful mainly for debugging purposes
     */
    std::string format(const ArgVector& av)
    {
      std::stringstream os;
      os << "[";
      if (av.empty())
      {
        os << "]";
      }
      else
      {
        os << av[0]();
        for (std::size_t i = 1; i < av.size(); ++i)
        {
          os << ", " << av[i]();
        }
        os << "]";
      }
      return os.str();
    }

    /**
     * @brief Map ArgType from Option to string used as part of help description
     * for an Option
     */
    inline const std::string& typeString(const Option& opt)
    {
      static const std::string labels[] =
          {"STRING", "REAL", "INT", "BOOL", "UNSPECIFIED"};
      return labels[static_cast<int>(opt.type)];
    }

    /**
     * @brief Stream insertion operator to print an Option specification and
     * parsed argument(s)
     *
     * @note This is useful mainly for debugging purposes
     */
    std::ostream& operator<<(std::ostream&             os,
                             const OptionArgumentPair& arg)
    {
      os << "  name     : " << arg.option.name[0];
      if (!arg.option.name[1].empty())
      {
        os << ", " << arg.option.name[1];
      }
      os << '\n';
      os << "  help     : " << arg.option.help << '\n';
      os << "  required : " << std::boolalpha << arg.option.required << '\n';
      os << "  type     : " << typeString(arg.option) << '\n';
      os << "  flag     : " << std::boolalpha << arg.option.flag << '\n';
      os << "  defaults : " << format(arg.option.defaults) << '\n';
      os << "  nargs    : " << arg.option.nargs << '\n';
      os << "  values   : " << format(arg.values);
      return os;
    }

    /// Implementation of CliArgs
    std::ostream& operator<<(std::ostream& os, const CliArgsImpl& args)
    {
      for (auto&& arg : args.table_)
      {
        os << arg << "\n\n";
      }
      for (auto&& [key, arg] : args.name_map_)
      {
        os << key << " =>\n"
           << *arg << "\n\n";
      }
      return os;
    }

    void CliArgsImpl::mapName(const std::string& name, OptionArgumentPair* arg)
    {
      if (!name_map_.try_emplace(name, arg).second)
      {
        auto msg = "CliArgs: attempt to add duplicate option name: \""
                   + name + "\".";
        Log::error() << msg << std::endl;
        throw std::runtime_error(msg);
      }
    }

    bool notANumber(const std::string& raw)
    {
      double v;
      return (std::stringstream(raw) >> v).fail();
    }

    void CliArgsImpl::setupArg(OptionArgumentPair& arg)
    {
      const auto& name = arg.option.name;
      assert(name[0].starts_with("--"));
      mapName(name[0].substr(2), &arg);
      if (!name[1].empty())
      {
        assert(name[1].starts_with("-"));
        assert(!name[1].starts_with("--"));
        if (!notANumber(name[1]))
        {
          auto msg = "CliArgs: Cannot use number as option name: \""
                     + name[1] + "\".";
          Log::error() << msg << std::endl;
          throw std::runtime_error(msg);
        }
        mapName(name[1].substr(1), &arg);
      }

      if (arg.option.flag)
      {
        arg.option.defaults = false;
        arg.option.nargs    = 0;
      }
    }

    CliArgsImpl::CliArgsImpl(std::initializer_list<Option> opts)
    {
      // Store and map options
      table_.reserve(opts.size() + 1);
      for (auto& opt : opts)
      {
        auto& arg = table_.emplace_back(opt);
        setupArg(arg);
      }

      // add "help" entry on behalf of the user
      if (!name_map_.count("help") && !name_map_.count("h"))
      {
        auto& help = table_.emplace_back(
            Option{.name = {"--help", "-h"},
                   .help = "Print this help message",
                   .flag = true});
        setupArg(help);
      }
    }

    const OptionArgumentPair* CliArgsImpl::findByName(const std::string& name)
        const
    {
      auto ret = name_map_.find(name);
      if (ret == name_map_.end())
      {
        return nullptr;
      }
      return ret->second;
    }

    bool looksLikeAnOption(const std::string& raw)
    {
      return raw.starts_with("-") && notANumber(raw);
    }

    OptionArgumentPair* CliArgsImpl::findByRawName(const std::string& raw)
    {
      if (raw.starts_with("--"))
      {
        return findByName(raw.substr(2));
      }
      else if (raw.starts_with("-"))
      {
        return findByName(raw.substr(1));
      }
      return nullptr;
    }

    void CliArgsImpl::parseArgs(int argc, const char* argv[])
    {
#if defined(_WIN32)
      app_name_ = std::filesystem::path(argv[0]).filename().string();
#else
      app_name_ = std::filesystem::path(argv[0]).filename();
#endif
      bool status = true;

      // Current argument (may involve multiple tokens)
      OptionArgumentPair* arg{nullptr};
      for (int i = 1; i < argc; ++i)
      {
        auto token = std::string(argv[i]);
        if (looksLikeAnOption(token))
        {
          auto* new_arg = findByRawName(token);
          if (!new_arg)
          {
            Log::error() << "CliArgs: unrecognized option \"" << token
                         << "\"\n";
            status = false;
            arg    = nullptr;
            continue;
          }
          if (arg && arg->values.size() != arg->option.nargs)
          {
            // option not given correct number of values
            Log::error() << "CliArgs: option \"" << arg->option.name[0]
                         << "\" requires " << arg->option.nargs
                         << " arguments; " << arg->values.size()
                         << " given.\n";
            status = false;
          }
          arg = new_arg;
          if (arg->option.flag)
          {
            // handle current arg as flag and move on
            arg->values = true;
          }
        }
        else if (arg)
        {
          if (arg->option.flag)
          {
            Log::error() << "CliArgs: flag option \"" << arg->option.name[0]
                         << "\" was given an argument (\"" << token << "\")\n";
            status = false;
          }
          // set value from current token
          arg->values.vec_.emplace_back(token);
        }
        else
        {
          // TODO: positional arguments
          Log::error() << "CliArgs: unrecognized option \"" << token
                       << "\"\n";
          status = false;
          arg    = nullptr;
        }
      }

      // Check for final flag argument
      if (arg)
      {
        // If we still have an active arg, it's either a flag (which is OK) or
        // an option which expects a value that hasn't been given
        if (!arg->option.flag && arg->values.size() != arg->option.nargs)
        {
          // option not given correct number of values
          Log::error() << "CliArgs: option \"" << arg->option.name[0]
                       << "\" requires " << arg->option.nargs << " arguments; "
                       << arg->values.size() << " given.\n";
          status = false;
        }
      }

      // check for help request
      if ((*this)["help"].as<bool>())
      {
        printHelp();
        exit(EXIT_SUCCESS);
      }

      // check that all required options have values
      for (const auto& arg : table_)
      {
        if (arg.option.required && arg.values.size() != arg.option.nargs)
        {
          Log::error() << "CliArgs: no input given for required option \""
                       << arg.option.name[0] << "\"\n";
          status = false;
        }
      }

      // print help and quit if errors have occurred
      if (!status)
      {
        printHelp();
        throw std::runtime_error("CliArgs: requirements not met");
      }
    }

    namespace
    {
      // Items in this anonymous namespace were stolen (and simplified) from
      // Boost.ProgramOptions

      /// A help line can't be longer than 80 characters
      inline constexpr unsigned maxLineLength = 80;

      /// Column for help descriptions will be half the total line length
      inline constexpr unsigned minHelpLength = maxLineLength / 2;

      /**
       * @brief Format string with option name(s) and type (if specified)
       */
      std::string formatFirstColumn(const Option& opt)
      {
        std::stringstream ss;
        ss << "  ";
        if (!opt.name[1].empty())
        {
          ss << opt.name[1] << ", ";
        }
        ss << opt.name[0];
        if (opt.type != ArgType::Unspecified)
        {
          ss << " (" << typeString(opt) << ")";
        }
        return ss.str();
      }

      /**
       * @brief Determine how much width the first column needs
       *
       * The width of the help column takes priority, and in cases where the
       * first column option string exceeds the arg column width, the help
       * description will simply be moved down a line.
       */
      unsigned getArgColumnWidth(const std::vector<OptionArgumentPair>& args)
      {
        const unsigned startOfHelpColumn = maxLineLength - minHelpLength;

        // Find the maximum width of the option column
        unsigned width{23};
        for (auto&& arg : args)
        {
          auto c1 = formatFirstColumn(arg.option);
          width   = std::max(width, static_cast<unsigned>(c1.size()));
        }

        // If first column is longer than the start of the help column,
        // we'll go to a new line
        width = std::min(width, startOfHelpColumn - 1);

        // add an additional space to improve readability
        ++width;
        return width;
      }

      /**
       * @brief Format a paragraph to the output stream
       *
       * This function attempts to take care of line breaks well
       */
      void printParagraph(std::ostream&    os,
                          std::string_view par,
                          unsigned         indent,
                          unsigned         lineLength)
      {
        // Through remainder of this function, 'lineLength' will be the length
        // available for characters, not including indent.
        assert(indent < lineLength);
        lineLength -= indent;

        if (par.size() < lineLength)
        {
          os << par;
        }
        else
        {
          const auto parEnd = par.cend();
          for (auto lineBegin = par.cbegin(); lineBegin != parEnd;)
          {
            // If line starts with space, but second character
            // is not space, remove the leading space.
            // We don't remove double spaces because those
            // might be intentianal.
            if ((*lineBegin == ' ') && ((lineBegin + 1 < parEnd) && (*(lineBegin + 1) != ' ')))
            {
              ++lineBegin;
            }

            unsigned remaining =
                static_cast<unsigned>(std::distance(lineBegin, parEnd));
            auto lineEnd =
                lineBegin + ((remaining < lineLength) ? remaining : lineLength);

            // prevent chopped words
            // Is lineEnd between two non-space characters?
            if ((*(lineEnd - 1) != ' ') && ((lineEnd < parEnd) && (*lineEnd != ' ')))
            {
              // find last ' ' in the second half of the current paragraph
              // line
              auto line      = std::string_view(lineBegin, lineEnd);
              auto lastSpace = line.find_last_of(' ');

              if (lastSpace != std::string_view::npos)
              {
                // is lastSpace within the second half of the current line
                if ((line.size() - lastSpace) < (lineLength / 2))
                {
                  lineEnd = lineBegin + lastSpace;
                }
              }
            }

            // write line to stream
            os << std::string_view(lineBegin, lineEnd);

            // more lines to follow?
            if (lineEnd != parEnd)
            {
              os << '\n';
              leftPad(os, indent);
            }

            lineBegin = lineEnd;
          }
        }
      }

      /**
       * @brief Print the help description of an option as a paragraph
       */
      void printOptionHelp(std::ostream&      os,
                           const std::string& desc,
                           unsigned           firstColWidth,
                           unsigned           lineLength)
      {
        // we need to use one char less per line to work correctly if actual
        // console has longer lines
        assert(lineLength > 1);
        if (lineLength > 1)
        {
          --lineLength;
        }

        // lineLength must be larger than firstColWidth
        // this assert may fail due to user error or environment conditions!
        assert(lineLength > firstColWidth);

        unsigned          padSize = 0;
        std::stringstream ss(desc);
        for (std::string par; std::getline(ss, par);)
        {
          leftPad(os, padSize);
          printParagraph(os, par, firstColWidth, lineLength);
          os << '\n';
          padSize = firstColWidth;
        }
      }

      /**
       * @brief Print the option name(s) [and type] followed by the help
       * description [and defaults]
       */
      void printOption(std::ostream& os,
                       const Option& opt,
                       unsigned      firstColWidth)
      {
        auto c1 = formatFirstColumn(opt);
        os << c1;

        if (!opt.help.empty())
        {
          unsigned padSize = firstColWidth - static_cast<unsigned>(c1.size());
          if (c1.size() >= firstColWidth)
          {
            os.put('\n'); // first column is too long, new line for help
            padSize = firstColWidth;
          }
          leftPad(os, padSize);

          auto desc = opt.help;
          if (!opt.flag && !opt.defaults.empty())
          {
            desc += "\ndefault: " + format(opt.defaults);
          }
          printOptionHelp(os, desc, firstColWidth, maxLineLength);
        }
      }
    } // namespace

    /**
     * @brief Format the portion of the usage expression for a single Option
     */
    std::string formatOptionUsage(const Option& opt)
    {
      std::stringstream ss;
      if (!opt.required)
      {
        ss << "[";
      }
      if (!opt.name[1].empty())
      {
        ss << opt.name[1];
      }
      else
      {
        ss << opt.name[0];
      }

      auto placeholder = toUpper(opt.name[0].substr(2));
      for (std::size_t i = 0; i < opt.nargs; ++i)
      {
        ss << " " << placeholder;
      }

      if (!opt.required)
      {
        ss << "]";
      }
      return ss.str();
    }

    void CliArgsImpl::printUsage(std::ostream& os) const
    {
      os << "Usage: " << app_name_;
      for (auto&& arg : table_)
      {
        os << " " << formatOptionUsage(arg.option);
      }
      os << "\n\n";
    }

    void CliArgsImpl::printHelp(std::ostream& os) const
    {
      printUsage(os);

      auto width = getArgColumnWidth(table_);

      for (auto&& arg : table_)
      {
        printOption(os, arg.option, width);
        os << "\n";
      }
    }

    const ArgVector& CliArgsImpl::operator[](const std::string& name) const
    {
      auto arg = findByName(name);
      if (!arg)
      {
        auto msg = "CliArgs: unrecognized option requested: \"" + name + "\"";
        Log::error() << msg << '\n';
        throw std::runtime_error(msg);
      }
      if (arg->values.empty() && !arg->option.defaults.empty())
      {
        return arg->option.defaults;
      }
      return arg->values;
    }
  } // namespace Utilities
} // namespace GridKit
