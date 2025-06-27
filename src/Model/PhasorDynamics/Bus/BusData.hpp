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

      /// Indices of the variables able to be monitored on this component
      enum class MonitorableVariables : size_t
      {
        VR,
        VI,
        VM,
        VA,
        MAXIMUM,
      };

      /// Set indicating the variables being monitored
      std::bitset<static_cast<
          std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::MAXIMUM)>
          monitored_variables;
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
        bd.bus_type = BusData<RealT, IdxT>::BusType::DEFAULT;
      }
      else if (string_class == "infinite_bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
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
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VR));
          }
          else if (monitored == "Vi")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VI));
          }
          else if (monitored == "Vm")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VM));
          }
          else if (monitored == "Va")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VA));
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
