#pragma once

#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>

#include <magic_enum/magic_enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     */
    template <typename real_type,
              typename index_type,
              typename parameters_enum,
              typename buses_enum,
              typename signal_inputs_enum,
              typename signal_outputs_enum,
              typename monitorable_variables_enum>
      requires std::is_enum_v<parameters_enum>
               && std::is_enum_v<buses_enum>
               && std::is_enum_v<signal_inputs_enum>
               && std::is_enum_v<signal_outputs_enum>
               && std::is_enum_v<monitorable_variables_enum>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;

      /// Index type
      using IdxT = index_type;

      /// Enumeration over parameters.
      using Parameters = parameters_enum;

      /// Enumeration over buses.
      using Buses = buses_enum;

      /// Enumeration over signal inputs.
      using SignalInputs = signal_inputs_enum;

      /// Enumeration over signal outputs.
      using SignalOutputs = signal_outputs_enum;

      /// Enumeration over monitorable variables.
      using MonitorableVariables = monitorable_variables_enum;

      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, std::variant<bool, RealT, IdxT>> parameters;

      /// Mapping of terminal attachments to bus identifiers
      std::map<Buses, IdxT> buses;

      /// Mapping of signal inputs to signal identifiers
      std::map<SignalInputs, IdxT> signal_inputs;

      /// Mapping of signal outputs to signal identifiers
      std::map<SignalOutputs, IdxT> signal_outputs;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      /// Disambiguation string for this device
      std::string disambiguation_string;
    };

    /**
     * @brief Numeric parameter `key` of `data`
     *
     * @throws std::invalid_argument if the parameter is missing or Boolean
     */
    template <typename DataT>
    typename DataT::RealT realParameter(const DataT& data, typename DataT::Parameters key)
    {
      const std::string name  = data.disambiguation_string + ": parameter " + std::string(magic_enum::enum_name(key));
      const auto        entry = data.parameters.find(key);
      if (entry == data.parameters.end())
      {
        throw std::invalid_argument(name + " is required");
      }
      if (const auto* real_value = std::get_if<typename DataT::RealT>(&entry->second))
      {
        return *real_value;
      }
      if (const auto* integer_value = std::get_if<typename DataT::IdxT>(&entry->second))
      {
        return static_cast<typename DataT::RealT>(*integer_value);
      }
      throw std::invalid_argument(name + " must be numeric");
    }

    /**
     * @brief Numeric parameter `key` of `data`, or `fallback` if it is not set
     */
    template <typename DataT>
    typename DataT::RealT realParameter(const DataT& data, typename DataT::Parameters key, typename DataT::RealT fallback)
    {
      if (!data.parameters.contains(key))
      {
        return fallback;
      }
      return realParameter(data, key);
    }
  } // namespace PhasorDynamics
} // namespace GridKit
