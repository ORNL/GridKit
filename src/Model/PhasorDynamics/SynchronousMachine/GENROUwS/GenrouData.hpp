/**
 * @file GenrouData.hpp
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
     * @brief Contains modeling data for a Genrou generator model.
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct GenrouData
    {
      RealT p0{0.0};    ///< Initial active power
      RealT q0{0.0};    ///< Initial reactive power
      RealT H{0.0};     ///< Rotor inertia
      RealT D{0.0};     ///< Damping coefficient
      RealT Ra{0.0};    ///< Winding resistance
      RealT Tdop{0.0};  ///< Open circuit direct axis transient time
      RealT Tdopp{0.0}; ///< Open circuit direct axis sub-transient time
      RealT Tqop{0.0};  ///< Open circuit quadrature axis transient
      RealT Tqopp{0.0}; ///< Open circuit quadrature axis sub-transient time
      RealT Xd{0.0};    ///< Direct axis synchronous reactance
      RealT Xdp{0.0};   ///< Direct axis transient reactance
      RealT Xdpp{0.0};  ///< Direct axis sub-transient reactance
      RealT Xq{0.0};    ///< Quadrature axis synchronous reactance
      RealT Xqp{0.0};   ///< Quadrature axis transient reactance
      RealT Xqpp{0.0};  ///< Quadrature axis sub-transient reactance
      RealT Xl{0.0};    ///< Stator leakage reactance
      RealT S10{0.0};   ///< Saturation factor at 1.0 pu flux
      RealT S12{0.0};   ///< Saturation factor at 1.2 pu flux

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

    template <typename RealT, typename IdxT>
    void from_json(const json& j, GenrouData<RealT, IdxT>& gd)
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
        else if (raw_parameter.key() == "Tdop")
        {
          raw_parameter.value().get_to(gd.Tdop);
        }
        else if (raw_parameter.key() == "Tdopp")
        {
          raw_parameter.value().get_to(gd.Tdopp);
        }
        else if (raw_parameter.key() == "Tqopp")
        {
          raw_parameter.value().get_to(gd.Tqopp);
        }
        else if (raw_parameter.key() == "Tqop")
        {
          raw_parameter.value().get_to(gd.Tqop);
        }
        else if (raw_parameter.key() == "Xd")
        {
          raw_parameter.value().get_to(gd.Xd);
        }
        else if (raw_parameter.key() == "Xdp")
        {
          raw_parameter.value().get_to(gd.Xdp);
        }
        else if (raw_parameter.key() == "Xdpp")
        {
          raw_parameter.value().get_to(gd.Xdpp);
        }
        else if (raw_parameter.key() == "Xq")
        {
          raw_parameter.value().get_to(gd.Xq);
        }
        else if (raw_parameter.key() == "Xqp")
        {
          raw_parameter.value().get_to(gd.Xqp);
        }
        else if (raw_parameter.key() == "Xqpp")
        {
          raw_parameter.value().get_to(gd.Xqpp);
        }
        else if (raw_parameter.key() == "Xl")
        {
          raw_parameter.value().get_to(gd.Xl);
        }
        else if (raw_parameter.key() == "S10")
        {
          raw_parameter.value().get_to(gd.S10);
        }
        else if (raw_parameter.key() == "S12")
        {
          raw_parameter.value().get_to(gd.S12);
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
                GenrouData<RealT, IdxT>::MonitorableVariables::IR));
          }
          else if (monitored == "ii")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenrouData<RealT, IdxT>::MonitorableVariables::II));
          }
          else if (monitored == "p")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenrouData<RealT, IdxT>::MonitorableVariables::P));
          }
          else if (monitored == "q")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenrouData<RealT, IdxT>::MonitorableVariables::Q));
          }
          else if (monitored == "delta")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenrouData<RealT, IdxT>::MonitorableVariables::DELTA));
          }
          else if (monitored == "omega")
          {
            gd.monitored_variables.set(static_cast<size_t>(
                GenrouData<RealT, IdxT>::MonitorableVariables::OMEGA));
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
