/**
 * @file BusData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

#include <optional>
#include <set>
#include <string>
#include <type_traits>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Indices of the variables able to be monitored on this component
    enum class BusMonitorableVariables : size_t
    {
      Vr,
      Vi,
      Vm,
      Va
    };

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
      std::string name; ///< A name given to this bus

      RealT Vr0{1.0}; ///< Initial value for the real bus voltage
      RealT Vi0{0.0}; ///< Initial value for the imaginary bus voltage

      IdxT bus_id{0}; ///< The unique ID of the bus

      /// Enumeration over the kinds of bus this data structure can be for
      enum class BusType
      {
        INVALID,
        DEFAULT,
        SLACK,
      };

      BusType bus_type{BusType::INVALID}; ///< The kind of bus this data is for

      RealT                v_base{1.0}; ///< Voltage base in volts
      std::optional<RealT> freq_base;   ///< Override for the system-wide base frequency
      std::optional<RealT> va_base;     ///< Override for the system-wide power base

      /// Alias
      using MonitorableVariables = BusMonitorableVariables;

      /// Set indicating the variables being monitored
      std::set<MonitorableVariables> monitored_variables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
