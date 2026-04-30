/**
 * @file GensalData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for a GENSAL generator model.
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a Gensal generator model
    enum class GensalParameters
    {
      p0,       ///< Initial active power
      q0,       ///< Initial reactive power
      H,        ///< Rotor inertia
      D,        ///< Damping coefficient
      Ra,       ///< Winding resistance
      Tdop,     ///< Open circuit direct axis transient time
      Tdopp,    ///< Open circuit direct axis sub-transient time
      Tqopp,    ///< Open circuit quadrature axis sub-transient time
      Xd,       ///< Direct axis synchronous reactance
      Xdp,      ///< Direct axis transient reactance
      Xdpp,     ///< Direct axis sub-transient reactance
      Xq,       ///< Quadrature axis synchronous reactance
      Xl,       ///< Stator leakage reactance
      S10,      ///< Saturation factor at 1.0 pu flux
      S12,      ///< Saturation factor at 1.2 pu flux
      mva_base, ///< MVA base of the gensal model
    };

    /// Ports for a Gensal generator model
    enum class GensalPorts
    {
      bus,   ///< Unique ID of the connecting bus
      pmech, ///< Unique ID of the signal providing mechanical power
      speed, ///< Unique ID of the signal receiving speed deviation
      efd,   ///< Unique ID of the signal providing exciter field voltage
    };

    /// Variables able to be monitored for a Gensal generator model
    enum class GensalMonitorableVariables
    {
      ir,
      ii,
      p,
      q,
      delta,
      omega,
      speed
    };

    /**
     * @brief Contains modeling data for a Gensal generator model.
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename RealT, typename IdxT>
    struct GensalData : public ComponentData<RealT,
                                             IdxT,
                                             GensalParameters,
                                             GensalPorts,
                                             GensalMonitorableVariables>
    {
      GensalData() = default;

      using Parameters           = GensalParameters;
      using Ports                = GensalPorts;
      using MonitorableVariables = GensalMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
