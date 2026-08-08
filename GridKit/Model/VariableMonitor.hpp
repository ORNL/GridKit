#pragma once

#include <algorithm>
#include <array>
#include <charconv>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <ranges>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <variant>
#include <vector>

#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Model
  {
    template <typename scalar_type>
    class VariableMonitorController;

    namespace VariableMonitorDetail
    {
      template <typename RealT>
      std::string formatReal(RealT value)
      {
        std::array<char, 128> buffer{};
        constexpr auto        precision = std::numeric_limits<RealT>::digits10 + 1;

        auto [ptr, ec] = std::to_chars(buffer.data(),
                                       buffer.data() + buffer.size(),
                                       value,
                                       std::chars_format::scientific,
                                       precision);

        if (ec == std::errc{})
        {
          return std::string(buffer.data(), ptr);
        }

        std::ostringstream os;
        os.precision(precision);
        os << std::scientific << value;
        return os.str();
      }
    } // namespace VariableMonitorDetail

    /**
     * @enum VariableMonitorFormat
     * Available formats for monitor output
     */
    enum class VariableMonitorFormat
    {
      CSV,  ///< CSV format
      JSON, ///< JSON format
      YAML  ///< YAML format
    };

    /**
     * @brief Selects which monitored columns a sink receives
     *
     * Patterns are matched against the full column name, so a component with
     * several monitored variables can be split across sinks.
     */
    class SinkFilter
    {
    public:
      SinkFilter() = default;

      explicit SinkFilter(std::vector<std::string> patterns)
        : patterns_(std::move(patterns))
      {
      }

      /// A filter with no patterns selects everything
      bool operator()(std::string_view column) const
      {
        if (patterns_.empty())
        {
          return true;
        }
        return std::ranges::any_of(patterns_,
                                   [column](const std::string& pattern)
                                   { return match(pattern, column); });
      }

    private:
      /**
       * @brief Glob match supporting `*` (any run) and `?` (any one character)
       *
       * Backtracks to the last `*` on mismatch, so no recursion is needed.
       */
      static bool match(std::string_view pattern, std::string_view text)
      {
        std::size_t p = 0, t = 0, star = std::string_view::npos, mark = 0;
        while (t < text.size())
        {
          if (p < pattern.size()
              && (pattern[p] == '?' || pattern[p] == text[t]))
          {
            ++p;
            ++t;
          }
          else if (p < pattern.size() && pattern[p] == '*')
          {
            star = p++;
            mark = t;
          }
          else if (star != std::string_view::npos)
          {
            p = star + 1;
            t = ++mark;
          }
          else
          {
            return false;
          }
        }
        while (p < pattern.size() && pattern[p] == '*')
        {
          ++p;
        }
        return p == pattern.size();
      }

      std::vector<std::string> patterns_;
    };

    /// Shared so that format tags stay cheap to copy
    using SinkFilterPtr = std::shared_ptr<const SinkFilter>;

    /**
     * @brief Abstract class for managing output of monitored variables
     *
     * This class is used for both the high-level control monitor and the
     * individual component and bus monitors.
     */
    class VariableMonitorBase
    {
      template <typename>
      friend class VariableMonitorController;

    public:
      /// Type used for dispatch
      struct Csv
      {
        /// Delimiter for CSV line output
        std::string   delim{","};
        /// Columns this sink accepts (null selects everything)
        SinkFilterPtr filter{};
      };

      /// Type used for dispatch
      struct Json
      {
        /// Implementation detail used to prevent a comma before the first block
        mutable bool  after_first{false};
        /// Columns this sink accepts (null selects everything)
        SinkFilterPtr filter{};
      };

      /// Type used for dispatch
      struct Yaml
      {
        /// Columns this sink accepts (null selects everything)
        SinkFilterPtr filter{};
      };

      /// Short alias for local use
      using Format = VariableMonitorFormat;

      /**
       * @brief Defines information necessary to create a monitor sink
       */
      struct SinkSpec
      {
        /// Output format
        Format                   format;
        /// Output file name (empty for stdout)
        std::string              file_name{};
        /// Delimiter (used only with CSV format currently)
        std::string              delim{","};
        /// Column name globs; empty selects every monitored column
        std::vector<std::string> include{};
      };

      virtual ~VariableMonitorBase()
      {
      }

      /**
       * @brief Is there nothing to monitor?
       */
      virtual bool empty() const = 0;

    protected:
      ///@{
      /**
       * @brief Print items relevant to the start of a file
       */
      virtual void appendHeader(std::string&, Csv) const = 0;

      virtual void appendHeader(std::string&, Json) const
      {
      }

      virtual void appendHeader(std::string&, Yaml) const
      {
      }

      ///@}

      ///@{
      /**
       * @brief Print monitored variables at current state
       */
      virtual void append(std::string&, Csv) const  = 0;
      virtual void append(std::string&, Json) const = 0;
      virtual void append(std::string&, Yaml) const = 0;

      ///@}

      ///@{
      /**
       * @brief Print items relevant to the end of a file
       */
      virtual void appendFooter(std::string&, Csv) const
      {
      }

      virtual void appendFooter(std::string&, Json) const
      {
      }

      virtual void appendFooter(std::string&, Yaml) const
      {
      }

      ///@}
    };

    template <typename eval_type, template <typename, typename> typename model_data_type>
    class VariableMonitor;

  } // namespace Model
} // namespace GridKit
