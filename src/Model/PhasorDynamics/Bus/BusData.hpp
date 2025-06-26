/**
 * @file BusData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
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
    using json = nlohmann::json;

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
      /// Enumeration over the kinsd of bus this data structure
      /// can be for
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

      /// The unique ID of the bus
      IdxT bus_id{0};

      /// The kind of bus this data is for
      BusType bus_type{BusType::Default};

      /// Voltage base in volts
      RealT v_base{1.0};

      /// Override for the system-wide base frequency
      std::optional<RealT> freq_base;

      /// Override for the system-wide power base
      std::optional<RealT> va_base;

      /// Indices of the variables able to be monitored on this
      /// component
      enum class MonitorableVariables : size_t
      {
        Vr,
        Vi,
        Vm,
        Va,
        Maximum,
      };

      /// Set indicating the variables which are monitored
      std::bitset<static_cast<std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::Maximum)> monitored_variables;
    };

    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bd)
    {
      j.at("name").get_to(bd.name);

      if (j.contains("init"))
      {
        for (auto& raw_parameter : j.at("init").items())
        {
          if (raw_parameter.key() == "Vi")
          {
            raw_parameter.value().get_to(bd.Vi0);
          }
          else if (raw_parameter.key() == "Vr")
          {
            raw_parameter.value().get_to(bd.Vr0);
          }
          else
          {
            throw "Invalid initial parameter";
          }
        }
      }

      j.at("number").get_to(bd.bus_id);

      auto string_class = j.at("class").get<std::string>();
      if (string_class == "bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::Default;
      }
      else if (string_class == "infinite_bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::Slack;
      }
      else
      {
        throw "Invalid bus class";
      }

      j.at("v_base").get_to(bd.v_base);

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(bd.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(bd.va_base);
      }

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto monitored = raw_monitored_variable.get<std::string>();
          if (monitored == "Vr")
          {
            bd.monitored_variables.set(static_cast<size_t>(BusData<RealT, IdxT>::MonitorableVariables::Vr));
          }
          else if (monitored == "Vi")
          {
            bd.monitored_variables.set(static_cast<size_t>(BusData<RealT, IdxT>::MonitorableVariables::Vi));
          }
          else if (monitored == "Vm")
          {
            bd.monitored_variables.set(static_cast<size_t>(BusData<RealT, IdxT>::MonitorableVariables::Vm));
          }
          else if (monitored == "Va")
          {
            bd.monitored_variables.set(static_cast<size_t>(BusData<RealT, IdxT>::MonitorableVariables::Va));
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
