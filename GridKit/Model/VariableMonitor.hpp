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
      void formatReal(std::string& out, RealT value)
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
          out.append(buffer.data(), ptr);
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
     * @enum VariableMonitorLayout
     * Available table layouts for monitor output
     */
    enum class VariableMonitorLayout
    {
      WIDE, ///< One row per sample with one column per variable
      LONG  ///< One record per sample and variable
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
      /// Short alias for local use
      using Layout = VariableMonitorLayout;

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
        /// Table layout
        Layout      layout{Layout::WIDE};
      };

      virtual ~VariableMonitorBase()
      {
      }

      /**
       * @brief Is there nothing to monitor?
       */
      virtual bool empty() const = 0;

      /// Returns the current value of a monitored variable
      using GetterT = std::function<double()>;

      /**
       * @brief Monitored variable and the getter of its current value
       */
      struct Entry
      {
        /// Device class, such as `Genrou`
        std::string_view device_class;
        /// Device identifier within its class
        std::string_view id;
        /// Variable name
        std::string_view variable;
        /// Returns the current value
        GetterT          value;
      };

    protected:
      /// Append one entry per monitored variable
      virtual void appendEntries(std::vector<Entry>&) const = 0;
    };

    template <typename eval_type, template <typename, typename> typename model_data_type>
    class VariableMonitor;

  } // namespace Model
} // namespace GridKit
