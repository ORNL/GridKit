#include <algorithm>
#include <cassert>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <GridKit/Utilities/CliOptions/CliArgs.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Utilities
  {
    using Log = ::GridKit::Utilities::Logger;

    std::ostream& operator<<(std::ostream& os, const ArgVector& av)
    {
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
      return os;
    }

    inline const std::string&
    typeString(const Argument& arg)
    {
      static const std::string labelArray[] =
          {"STR", "REAL", "INT", "BOOL", "UNSPECIFIED"};
      return labelArray[static_cast<int>(arg.type)];
    }

    std::ostream& operator<<(std::ostream& os, const Argument& arg)
    {
      os << "name: " << arg.name[0];
      if (arg.name.size() == 2)
      {
        os << ", " << arg.name[1] << '\n';
      }
      os << "help: " << arg.help << '\n';
      os << "required: " << std::boolalpha << arg.required << '\n';
      os << "type: " << typeString(arg) << '\n';
      os << "flag: " << std::boolalpha << arg.flag << '\n';
      os << "defaults: " << arg.defaults << '\n';
      os << "nargs: " << arg.nargs << '\n';
      os << "values: " << arg.values << '\n';
      return os;
    }

    std::ostream& operator<<(std::ostream& os, const CliArgs& args)
    {
      for (auto&& arg : args.table_)
      {
        os << arg << '\n';
      }
      for (auto&& [key, arg] : args.name_map_)
      {
        os << key << " =>\n"
           << *arg << '\n';
      }
      return os;
    }

    void CliArgs::mapName(const std::string& name, Argument* arg)
    {
      if (!name_map_.try_emplace(name, arg).second)
      {
        Log::error() << "CliArgs: attempt to add duplicate option name: \""
                     << name << "\".\n";
      }
    }

    void CliArgs::setupArg(Argument& arg)
    {
      assert(arg.name.size() == 1 || arg.name.size() == 2);
      assert(arg.name.at(0).starts_with("--"));
      mapName(arg.name[0].substr(2), &arg);
      if (arg.name.size() == 2)
      {
        assert(arg.name[1].starts_with("-"));
        assert(!arg.name[1].starts_with("--"));
        mapName(arg.name[1].substr(1), &arg);
      }

      if (arg.flag)
      {
        arg.values.vec.emplace_back(false);
        arg.nargs = 0;
      }
    }

    CliArgs::CliArgs(std::initializer_list<Argument> args)
      : table_(args)
    {
      table_.reserve(table_.size() + 1);
      for (auto&& arg : table_)
      {
        setupArg(arg);
      }
      if (!name_map_.count("help") && !name_map_.count("h"))
      {
        auto& help = table_.emplace_back(
            Argument{.name = {"--help", "-h"},
                     .help = "Print this help message",
                     .flag = true});
        setupArg(help);
      }
    }

    Argument* CliArgs::findArg(const std::string& raw)
    {
      if (raw.starts_with("--"))
      {
        return name_map_.find(raw.substr(2))->second;
      }
      else if (raw.starts_with("-"))
      {
        return name_map_.find(raw.substr(1))->second;
      }
      return nullptr;
    }

    void CliArgs::parseArgs(int argc, const char* argv[])
    {
      app_name_ = std::filesystem::path(argv[0]).filename();

      Argument* arg{nullptr};
      for (int i = 1; i < argc; ++i)
      {
        auto token = std::string(argv[i]);
        if (!arg)
        {
          arg = findArg(token);
          if (!arg)
          {
            // TODO: positional argument(s)?
            Log::error() << "CliArgs: unrecognized option \"" << token
                         << "\"\n";
          }
        }
        else
        {
          auto* new_arg = findArg(token);
          if (new_arg)
          {
            if (arg->flag)
            {
              // handle current arg as flag and move on
              arg->values.vec[0] = true;
              arg                = new_arg;
            }
            else if (arg->values.vec.size() < arg->nargs)
            {
              // option not given correct number of values
              Log::error() << "CliArgs: option \"" << arg->name[0]
                           << "\" requires " << arg->nargs << " arguments; "
                           << arg->values.vec.size() << " given.\n";
            }
          }
          else
          {
            // set value from current token
            arg->values.vec.emplace_back(token);
            if (arg->values.vec.size() == arg->nargs)
            {
              // reset for next arg
              arg = nullptr;
            }
          }
        }
      }
      if (arg)
      {
        if (arg->flag)
        {
          arg->values.vec[0] = true;
        }
        else
        {
          // option not given correct number of values
          Log::error() << "CliArgs: option \"" << arg->name[0]
                       << "\" requires " << arg->nargs << " arguments; "
                       << arg->values.vec.size() << " given.\n";
        }
      }

      // check for help request
      if ((*this)["help"]().as<bool>())
      {
        printHelp();
        exit(EXIT_SUCCESS);
      }

      // check that all required arguments have values
      for (const auto& arg : table_)
      {
        bool err = false;
        if (arg.required && arg.values.vec.size() != arg.nargs)
        {
          Log::error() << "CliArgs: no input given for required argument \""
                       << arg.name[0] << "\"\n";
          err = true;
        }
        if (err)
        {
          printHelp();
          throw std::runtime_error("CliArgs: requirements not met");
        }
      }
    }

    namespace
    {
      // This is stolen (and simplified) from Boost.ProgramOptions
      inline constexpr unsigned defaultLineLength = 80;

      unsigned
      getLineLength()
      {
        return defaultLineLength;
      }

      unsigned
      getMinHelpLength()
      {
        return defaultLineLength / 2;
      }

      std::string
      formatFirstColumn(const Argument& arg)
      {
        std::stringstream ss;
        ss << "  ";
        if (arg.name.size() == 2)
        {
          ss << arg.name[1] << ", ";
        }
        ss << arg.name[0];
        if (arg.type != Argument::Type::Unspecified)
        {
          ss << " (" << typeString(arg) << ")";
        }
        return ss.str();
      }

      unsigned
      getArgColumnWidth(
          const std::vector<Argument>& args,
          const unsigned               lineLength,
          const unsigned               minHelpLength)
      {
        const unsigned startOfHelpColumn = lineLength - minHelpLength;

        // Find the maximum width of the option column
        unsigned width{23};
        for (auto&& arg : args)
        {
          auto c1 = formatFirstColumn(arg);
          width   = std::max(width, static_cast<unsigned>(c1.size()));
        }

        // If first column is longer than the start of the help column,
        // we'll go to a new line
        width = std::min(width, startOfHelpColumn - 1);

        // add an additional space to improve readability
        ++width;
        return width;
      }

      void
      leftPad(std::ostream& os, unsigned width)
      {
        os << std::setfill(' ') << std::setw(width) << "";
      }

      void
      printParagraph(std::ostream& os, std::string_view par, unsigned indent, unsigned lineLength)
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

      void
      printArgHelp(std::ostream& os, const std::string& desc, unsigned firstColWidth, unsigned lineLength)
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

      void
      printArg(std::ostream& os, const Argument& arg, unsigned firstColWidth, unsigned lineLength)
      {
        auto c1 = formatFirstColumn(arg);
        os << c1;

        if (!arg.help.empty())
        {
          unsigned padSize = firstColWidth - static_cast<unsigned>(c1.size());
          if (c1.size() >= firstColWidth)
          {
            os.put('\n'); // first column is too long, new line for help
            padSize = firstColWidth;
          }
          leftPad(os, padSize);

          auto desc = arg.help;
          if (!arg.defaults.empty())
          {
            desc += "\ndefault: " + (std::stringstream() << arg.defaults).str();
          }
          printArgHelp(os, desc, firstColWidth, lineLength);
        }
      }
    } // namespace

    std::string toUpper(std::string str)
    {
      std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c)
                     { return std::toupper(c); });
      return str;
    }

    std::string formatArgUsage(const Argument& arg)
    {
      std::stringstream ss;
      if (!arg.required)
      {
        ss << "[";
      }
      if (arg.name.size() == 2)
      {
        ss << arg.name[1];
      }
      else
      {
        ss << arg.name[0];
      }

      auto placeholder = toUpper(arg.name[0].substr(2));
      for (std::size_t i = 0; i < arg.nargs; ++i)
      {
        ss << " " << placeholder;
      }

      if (!arg.required)
      {
        ss << "]";
      }
      return ss.str();
    }

    void CliArgs::printHelp(std::ostream& os) const
    {
      os << "Usage: " << app_name_;
      for (auto&& arg : table_)
      {
        os << " " << formatArgUsage(arg);
      }
      os << "\n\n";

      auto lineLength    = getLineLength();
      auto minHelpLength = getMinHelpLength();
      auto width         = getArgColumnWidth(table_, lineLength, minHelpLength);

      for (auto&& arg : table_)
      {
        printArg(os, arg, width, lineLength);
        os << "\n";
      }
    }

  } // namespace Utilities
} // namespace GridKit
