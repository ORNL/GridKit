/**
 * @file BranchData.hpp
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
    using json = nlohmann::json;

    /**
     * @brief Contains modeling data for a Branch
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct BranchData
    {
      /// Line series resistance
      RealT R{0.0};

      /// Line series reactance
      RealT X{0.0};

      /// Line shunt conductance
      RealT G{0.0};

      /// Line shunt charging
      RealT B{0.0};

      /// Unique ID of bus 1
      IdxT bus1_id{0};

      /// Unique ID of bus 2
      IdxT bus2_id{0};

      /// Override for the system-wide base frequency
      std::optional<RealT> freq_base;

      /// Override for the system-wide power base
      std::optional<RealT> va_base;

      /// Disambiguation string for this device
      std::string disambiguation_string;

      /// Indices of the variables able to be monitored in the bitset
      enum class MonitorableVariables : size_t
      {
        Ir1,
        Ii1,
        Im1,
        P1,
        Q1,
        Ir2,
        Ii2,
        Im2,
        P2,
        Q2,
        Maximum,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::Maximum)> monitored_variables;
    };

    template <typename RealT, typename IdxT>
    void from_json(const json& j, BranchData<RealT, IdxT>& bd)
    {
      for (auto& raw_parameter : j.at("params").items())
      {
        if (raw_parameter.key() == "R")
        {
          raw_parameter.value().get_to(bd.R);
        }
        else if (raw_parameter.key() == "X")
        {
          raw_parameter.value().get_to(bd.X);
        }
        else if (raw_parameter.key() == "G")
        {
          raw_parameter.value().get_to(bd.G);
        }
        else if (raw_parameter.key() == "B")
        {
          raw_parameter.value().get_to(bd.B);
        }
        else
        {
          throw "Invalid initial parameter";
        }
      }

      for (auto& raw_port : j.at("ports").items())
      {
        if (raw_port.key() == "bus1")
        {
          raw_port.value().get_to(bd.bus1_id);
        }
        else if (raw_port.key() == "bus2")
        {
          raw_port.value().get_to(bd.bus2_id);
        }
        else
        {
          throw "Invalid port mapping";
        }
      }

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(bd.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(bd.va_base);
      }

      j.at("id").get_to(bd.disambiguation_string);

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto monitored = raw_monitored_variable.get<std::string>();
          if (monitored == "ir1")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Ir1));
          }
          else if (monitored == "ii1")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Ii1));
          }
          else if (monitored == "im1")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Im1));
          }
          else if (monitored == "p1")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::P1));
          }
          else if (monitored == "q1")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Q1));
          }
          else if (monitored == "ir2")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Ir2));
          }
          else if (monitored == "ii2")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Ii2));
          }
          else if (monitored == "im2")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Im2));
          }
          else if (monitored == "p2")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::P2));
          }
          else if (monitored == "q2")
          {
            bd.monitored_variables.set(static_cast<size_t>(BranchData<RealT, IdxT>::MonitorableVariables::Q2));
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
