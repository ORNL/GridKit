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
      void appendCsvReal(std::string& out, RealT value)
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
      virtual void printHeader(std::ostream&, Csv) const = 0;

      virtual void printHeader(std::ostream&, Json) const
      {
      }

      virtual void printHeader(std::ostream&, Yaml) const
      {
      }

      virtual void appendCsvHeader(std::string& out, Csv csv) const
      {
        std::ostringstream os;
        printHeader(os, csv);
        out += os.str();
      }

      ///@}

      ///@{
      /**
       * @brief Print monitored variables at current state
       */
      virtual void print(std::ostream&, Csv) const  = 0;
      virtual void print(std::ostream&, Json) const = 0;
      virtual void print(std::ostream&, Yaml) const = 0;

      virtual void appendCsv(std::string& out, Csv csv) const
      {
        std::ostringstream os;
        constexpr auto     max_prec = std::numeric_limits<double>::digits10 + 1;
        os.precision(max_prec);
        os << std::scientific;
        print(os, csv);
        out += os.str();
      }

      ///@}

      ///@{
      /**
       * @brief Print items relevant to the end of a file
       */
      virtual void printFooter(std::ostream&, Csv) const
      {
      }

      virtual void printFooter(std::ostream&, Json) const
      {
      }

      virtual void printFooter(std::ostream&, Yaml) const
      {
      }

      ///@}
    };

    template <typename EvalT, template <typename, typename> typename DataT>
    class VariableMonitor;

  } // namespace Model
} // namespace GridKit
