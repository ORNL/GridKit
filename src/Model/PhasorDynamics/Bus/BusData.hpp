/**
 * @file BusData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

#include <bitset>
#include <optional>
#include <type_traits>
#include <string>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a Bus
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct BusData
    {
      /// Enumeration over the kinds of bus this data structure can be for
      enum class BusType
      {
        Default,
        Slack,
      };

      /// A name given to this bus
      std::string name;

      /// Initial value for the real bus voltage
      RealT Vr0{1.0};

      /// Initial value for the imaginary bus voltage
      RealT Vi0{0.0};

      /// The unique ID of bus 1
      IdxT bus_id{0};

      /// The kind of bus this data is for
      BusType bus_type{BusType::Default};

      /// Voltage base in volts
      RealT v_base{1.0};

      /// Override for the system-wide base frequency
      std::optional<RealT> freq_base;

      /// Override for the system-wide power base
      std::optional<RealT> va_base;

      /// Indices of the variables able to be monitored in the bitset
      enum class MonitorableVariables : size_t
      {
        Vr,
        Vi,
        Vm,
        Va,
        Maximum,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::Maximum) - 1> monitored_variables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
