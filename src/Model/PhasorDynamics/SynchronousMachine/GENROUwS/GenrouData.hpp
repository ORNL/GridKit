/**
 * @file GenrouData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
 *
 */
#pragma once

#include <Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a Genrou generator model
    enum class GenrouParameters
    {
      /// Initial active power
      p0,

      /// Initial reactive power
      q0,

      /// Rotor inertia
      H,

      /// Damping coefficient
      D,

      /// Winding resistance
      Ra,

      /// Open circuit direct axis transient time
      Tdop,

      /// Open circuit direct axis sub-transient time
      Tdopp,

      /// Open circuit quadrature axis transient
      Tqop,

      /// Open circuit quadrature axis sub-transient time
      Tqopp,

      /// Direct axis synchronous reactance
      Xd,

      /// Direct axis transient reactance
      Xdp,

      /// Direct axis sub-transient reactance
      Xdpp,

      /// Quadrature axis synchronous reactance
      Xq,

      /// Quadrature axis transient reactance
      Xqp,

      /// Quadrature axis sub-transient reactance
      Xqpp,

      /// Stator leakage reactance
      Xl,

      /// Saturation factor at 1.0 pu flux
      S10,

      /// Saturation factor at 1.2 pu flux
      S12,
    };

    /// Ports for a Genrou generator model
    enum class GenrouPorts
    {
      /// Unique ID of the connecting bus
      bus,

      /// Unique ID of the bus providing the exciter signal
      exciter_signal,

      /// Unique ID of the bus providing the governor signal
      governor_signal,
    };

    /// Variables able to be monitored for a Genrou generator model
    enum class GenrouMonitorableVariables
    {
      ir,
      ii,
      p,
      q,
      delta,
      omega,
    };

    /**
     * @brief Contains modeling data for a Genrou generator model.
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename RealT, typename IdxT>
    struct GenrouData : public ComponentData<RealT,
                                             IdxT,
                                             GenrouParameters,
                                             GenrouPorts,
                                             GenrouMonitorableVariables>
    {
      GenrouData() = default;

      using Parameters           = GenrouParameters;
      using Ports                = GenrouPorts;
      using MonitorableVariables = GenrouMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
