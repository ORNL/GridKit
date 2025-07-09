/**
 * @file GenClassicalData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
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
     * @brief Contains modeling data for a GenClassical generator model.
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct GenClassicalData
    {
      IdxT  unit_id{0}; ///< Unique component ID
      RealT p0{0.0};    ///< Initial active power
      RealT q0{0.0};    ///< Initial reactive power
      RealT H{0.0};     ///< Rotor inertia
      RealT D{0.0};     ///< Damping coefficient
      RealT Ra{0.0};    ///< Winding resistance
      RealT Xdp{0.0};   ///< Direct axis transient reactance

      IdxT                bus_id{0};       ///< Unique ID of the connecting bus
      std::optional<IdxT> exciter_signal;  ///< Unique ID of the bus providing the exciter signal
      std::optional<IdxT> governor_signal; ///< Unique ID of the bus providing the governor signal

      std::optional<RealT> freq_base; ///< Override for the system-wide base frequency
      std::optional<RealT> va_base;   ///< Override for the system-wide power base

      std::string disambiguation_string; ///< Disambiguation string for this device

      /// Indices of the variables able to be monitored in the bitset
      enum class MonitorableVariables : size_t
      {
        IR,
        II,
        P,
        Q,
        DELTA,
        OMEGA,
        MAXIMUM,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<
          std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::MAXIMUM)>
          monitored_variables;
    };

    /// JSON parser function implementation for the `GenClassicalData` type
    ///
    /// See the `README.md` in `src/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, GenClassicalData<RealT, IdxT>& gd)
    {
      for (auto& raw_parameter : j.at("params").items())
      {
        if (raw_parameter.key() == "p0")
        {
          raw_parameter.value().get_to(gd.p0);
        }
        else if (raw_parameter.key() == "q0")
        {
          raw_parameter.value().get_to(gd.q0);
        }
        else if (raw_parameter.key() == "H")
        {
          raw_parameter.value().get_to(gd.H);
        }
        else if (raw_parameter.key() == "D")
        {
          raw_parameter.value().get_to(gd.D);
        }
        else if (raw_parameter.key() == "Ra")
        {
          raw_parameter.value().get_to(gd.Ra);
        }
        else if (raw_parameter.key() == "Xdp")
        {
          raw_parameter.value().get_to(gd.Xdp);
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
          raw_port.value().get_to(gd.bus_id);
        }
        else if (raw_port.key() == "governor_signal")
        {
          raw_port.value().get_to(gd.governor_signal);
        }
        else if (raw_port.key() == "exciter_signal")
        {
          raw_port.value().get_to(gd.exciter_signal);
        }
        else
        {
          throw "Invalid port mapping";
        }
      }

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(gd.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(gd.va_base);
      }

      j.at("id").get_to(gd.disambiguation_string);

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto monitored = raw_monitored_variable.get<std::string>();
          if (monitored == "ir")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenClassicalData<RealT, IdxT>::MonitorableVariables::IR));
          }
          else if (monitored == "ii")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenClassicalData<RealT, IdxT>::MonitorableVariables::II));
          }
          else if (monitored == "p")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenClassicalData<RealT, IdxT>::MonitorableVariables::P));
          }
          else if (monitored == "q")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenClassicalData<RealT, IdxT>::MonitorableVariables::Q));
          }
          else if (monitored == "delta")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenClassicalData<RealT, IdxT>::MonitorableVariables::DELTA));
          }
          else if (monitored == "omega")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenClassicalData<RealT, IdxT>::MonitorableVariables::OMEGA));
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
