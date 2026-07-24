#pragma once

#include <array>
#include <string>

namespace GridKit
{
  namespace Testing
  {

    /**
     * @brief Conveniently create a command-line string
     */
    template <int NArgs>
    struct CommandLine
    {
      CommandLine() = default;

      CommandLine(const std::array<std::string, NArgs>& args)
        : args_(args)
      {
        for (std::size_t i = 0; i < static_cast<std::size_t>(argc); ++i)
        {
          argv[i] = args_[i].c_str();
        }
      }

      template <typename... TT>
      CommandLine(const std::string& arg0, TT&&... args)
        : CommandLine(std::array<std::string, NArgs>{arg0, args...})
      {
      }

      int         argc{NArgs};
      const char* argv[NArgs + 1]{nullptr};

    private:
      const std::array<std::string, NArgs> args_;
    };

    template <typename... TT>
    CommandLine(const std::string&, TT...) -> CommandLine<1 + sizeof...(TT)>;

  } // namespace Testing
} // namespace GridKit
