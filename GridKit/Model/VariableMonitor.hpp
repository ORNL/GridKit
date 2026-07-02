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
      CSV,         ///< CSV format
      JSON,        ///< JSON format
      YAML,        ///< YAML format
      ARROW,       ///< Apache Arrow IPC file format (Feather v2)
      ARROW_STREAM ///< Apache Arrow IPC stream format
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

      /// Type used for dispatch (Arrow IPC output; carries no state)
      struct Arrow
      {
      };

      /// Short alias for local use
      using Format = VariableMonitorFormat;

      /**
       * @brief Defines information necessary to create a monitor sink
       */
      struct SinkSpec
      {
        /// Output format
        Format      format;
        /// Output file name (empty for stdout)
        std::string file_name{};
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
      virtual void appendHeader(std::string&, Csv) const = 0;

      virtual void appendHeader(std::string&, Json) const
      {
      }

      virtual void appendHeader(std::string&, Yaml) const
      {
      }

      /**
       * @brief Collect column names for Arrow output (one per monitored
       * variable)
       */
      virtual void appendHeader(std::vector<std::string>&, Arrow) const = 0;

      ///@}

      ///@{
      /**
       * @brief Print monitored variables at current state
       */
      virtual void append(std::string&, Csv) const  = 0;
      virtual void append(std::string&, Json) const = 0;
      virtual void append(std::string&, Yaml) const = 0;

      /**
       * @brief Collect current values for Arrow output; order matches the
       * column names
       */
      virtual void append(std::vector<double>&, Arrow) const = 0;

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
