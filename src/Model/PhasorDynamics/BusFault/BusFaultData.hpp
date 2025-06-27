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

#include <nlohmann/json.hpp>

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
      RealT R{0.0};        ///< Short to ground resistance
      RealT X{0.0};        ///< Short to ground reactance
      bool  status{false}; ///< If the fault has happened

      std::optional<IdxT> control_signal; ///< Unique ID of the bus providing a control signal
      IdxT                bus_id{0};      ///< Unique ID of the bus where the fault occurs

      std::optional<RealT> freq_base; ///< Override for the system-wide base frequency
      std::optional<RealT> va_base;   ///< Override for the system-wide power base

      std::string disambiguation_string; ///< Disambiguation string for this device

      /// Indices of the variables able to be monitored in the bitset
      enum class MonitorableVariables : size_t
      {
        STATE,
        IR,
        II,
        MAXIMUM,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<
          std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::MAXIMUM)>
          monitored_variables;
    };

    /// JSON parser function implementation for the `BusFaultData` type
    ///
    /// See the `README.md` in `src/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusFaultData<RealT, IdxT>& bf)
    {
      for (auto& raw_parameter : j.at("params").items())
      {
        if (raw_parameter.key() == "R")
        {
          raw_parameter.value().get_to(bf.R);
        }
        else if (raw_parameter.key() == "X")
        {
          raw_parameter.value().get_to(bf.X);
        }
        else if (raw_parameter.key() == "state0")
        {
          raw_parameter.value().get_to(bf.status);
        }
        else
        {
          throw "Invalid initial parameter";
        }
      }

      for (auto& raw_port : j.at("ports").items())
      {
        if (raw_port.key() == "bus")
        {
          raw_port.value().get_to(bf.bus_id);
        }
        else if (raw_port.key() == "control_signal")
        {
          raw_port.value().get_to(bf.control_signal);
        }
        else
        {
          throw "Invalid port mapping";
        }
      }

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(bf.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(bf.va_base);
      }

      j.at("id").get_to(bf.disambiguation_string);

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto monitored = raw_monitored_variable.get<std::string>();
          if (monitored == "state")
          {
            bf.monitored_variables.set(static_cast<size_t>(
                BusFaultData<RealT, IdxT>::MonitorableVariables::STATE));
          }
          else if (monitored == "ir")
          {
            bf.monitored_variables.set(static_cast<size_t>(
                BusFaultData<RealT, IdxT>::MonitorableVariables::IR));
          }
          else if (monitored == "ii")
          {
            bf.monitored_variables.set(static_cast<size_t>(
                BusFaultData<RealT, IdxT>::MonitorableVariables::II));
          }
          else
          {
            throw "Invalid monitored variable";
          }
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
