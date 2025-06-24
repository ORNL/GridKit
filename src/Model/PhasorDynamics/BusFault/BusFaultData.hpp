/**
 * @file BusFaultData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for short-to-ground fault
 *
 */
#pragma once

#include <bitset>
#include <optional>
#include <string>
#include <type_traits>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a short-to-ground fault
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct BusFaultData
    {
      /// Short to ground resistance
      RealT R{0.0};

      /// Short to ground reactance
      RealT X{0.0};

      /// If the fault has happened
      bool status{false};

      /// Unique ID of the bus where the fault occurs
      IdxT bus_id{0};

      /// Override for the system-wide base frequency
      std::optional<RealT> freq_base;

      /// Override for the system-wide power base
      std::optional<RealT> va_base;

      /// Disambiguation string for this device
      std::string disambiguation_string;

      /// Indices of the variables able to be monitored in the bitset
      enum class MonitorableVariables : size_t
      {
        State,
        Ir,
        Ii,
        Maximum,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::Maximum)> monitored_variables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
