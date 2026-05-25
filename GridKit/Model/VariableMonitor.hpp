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
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <type_traits>
#include <variant>
#include <vector>

#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Model
  {
    template <typename ScalarT>
    class VariableMonitorController;

    namespace VariableMonitorDetail
    {
      template <typename RealT>
      void appendReal(std::string& out, RealT value)
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
          out.append(buffer.data(), static_cast<std::size_t>(ptr - buffer.data()));
          return;
        }

        std::ostringstream os;
        os.precision(precision);
        os << std::scientific << value;
        out += os.str();
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
        std::string delim{","};
      };

      /// Type used for dispatch
      struct Json
      {
        /// Implementation detail used to prevent a comma before the first block
        mutable bool after_first{false};
      };

      /// Type used for dispatch
      struct Yaml
      {
      };

      /// Short alias for local use
      using Format = VariableMonitorFormat;

      /**
       * @brief Defines information necessary to create a monitor sink
       */
      struct SinkSpec
      {
        /// Output file name (empty for stdout)
        std::string file_name;
        /// Output format
        Format      format;
        /// Delimiter (used only with CSV format currently)
        std::string delim{","};
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
      virtual void printHeader(std::ostream& os, Csv csv) const
      {
        std::string out;
        appendHeader(out, csv);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void printHeader(std::ostream& os, Json json) const
      {
        std::string out;
        appendHeader(out, json);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void printHeader(std::ostream& os, Yaml yaml) const
      {
        std::string out;
        appendHeader(out, yaml);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

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
      virtual void print(std::ostream& os, Csv csv) const
      {
        std::string out;
        append(out, csv);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void print(std::ostream& os, Json json) const
      {
        std::string out;
        append(out, json);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void print(std::ostream& os, Yaml yaml) const
      {
        std::string out;
        append(out, yaml);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void append(std::string&, Csv) const  = 0;
      virtual void append(std::string&, Json) const = 0;
      virtual void append(std::string&, Yaml) const = 0;

      ///@}

      ///@{
      /**
       * @brief Print items relevant to the end of a file
       */
      virtual void printFooter(std::ostream& os, Csv csv) const
      {
        std::string out;
        appendFooter(out, csv);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void printFooter(std::ostream& os, Json json) const
      {
        std::string out;
        appendFooter(out, json);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

      virtual void printFooter(std::ostream& os, Yaml yaml) const
      {
        std::string out;
        appendFooter(out, yaml);
        os.write(out.data(), static_cast<std::streamsize>(out.size()));
      }

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

    template <typename EvalT, template <typename, typename> typename DataT>
    class VariableMonitor;

  } // namespace Model
} // namespace GridKit
