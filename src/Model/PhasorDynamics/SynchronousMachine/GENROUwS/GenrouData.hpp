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
      /// Unique unit ID
      IdxT unit_id{0};

      /// Initial active power
      RealT p0{0.0};

      /// Initial reactive power
      RealT q0{0.0};

      /// Rotor inertia
      RealT H{0.0};

      /// Damping coefficient
      RealT D{0.0};

      /// Winding resistance
      RealT Ra{0.0};

      /// Open circuit direct axis transient time
      RealT Tdop{0.0};

      /// Open circuit direct axis sub-transient time
      RealT Tdopp{0.0};

      /// Open circuit quadrature axis transient
      RealT Tqop{0.0};

      /// Open circuit quadrature axis sub-transient time
      RealT Tqopp{0.0};

      /// Direct axis synchronous reactance
      RealT Xd{0.0};

      /// Direct axis transient reactance
      RealT Xdp{0.0};

      /// Direct axis sub-transient reactance
      RealT Xdpp{0.0};

      /// Quadrature axis synchronous reactance
      RealT Xq{0.0};

      /// Quadrature axis transient reactance
      RealT Xqp{0.0};

      /// Quadrature axis sub-transient reactance
      RealT Xqpp{0.0};

      /// Stator leakage reactance
      RealT Xl{0.0};

      /// Saturation factor at 1.0 pu flux
      RealT S10{0.0};

      /// Saturation factor at 1.2 pu flux
      RealT S12{0.0};

      /// Unique ID of the connecting bus
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
        Ir,
        Ii,
        P,
        Q,
        Delta,
        Omega,
        Maximum,
      };

      /// Set of variables being monitored
      std::bitset<static_cast<std::underlying_type_t<MonitorableVariables>>(MonitorableVariables::Maximum) - 1> monitored_variables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
