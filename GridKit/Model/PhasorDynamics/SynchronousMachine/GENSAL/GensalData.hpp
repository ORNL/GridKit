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
      p0,    ///< Legacy fallback for the p input
      q0,    ///< Legacy fallback for the q input
      H,     ///< Rotor inertia
      D,     ///< Damping coefficient
      Ra,    ///< Winding resistance
      Tdop,  ///< Open circuit direct axis transient time
      Tdopp, ///< Open circuit direct axis sub-transient time
      Tqopp, ///< Open circuit quadrature axis sub-transient time
      Xd,    ///< Direct axis synchronous reactance
      Xdp,   ///< Direct axis transient reactance
      Xdpp,  ///< Direct axis sub-transient reactance
      Xq,    ///< Quadrature axis synchronous reactance
      Xl,    ///< Stator leakage reactance
      S10,   ///< Saturation factor at 1.0 pu flux
      S12,   ///< Saturation factor at 1.2 pu flux
      mva,   ///< MVA base of the gensal model
    };

    /// Buses for a Gensal generator model
    enum class GensalBuses : size_t
    {
      bus, ///< Unique ID of the connecting bus
      SIZE
    };

    /// Signal inputs for a Gensal generator model
    enum class GensalSignalInputs : size_t
    {
      pmech,  ///< Unique ID of the signal providing mechanical power
      efd,    ///< Unique ID of the signal providing exciter field voltage
      p,      ///< Initial active-power injection
      q,      ///< Initial reactive-power injection
      online, ///< In-service status (zero is offline, nonzero is online)
      SIZE
    };

    /// Signal outputs for a Gensal generator model
    enum class GensalSignalOutputs : size_t
    {
      speed, ///< Unique ID of the signal receiving speed deviation
      SIZE
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
      speed,
      Eqp,
      psidp,
      psiqpp,
      psidpp,
      vd,
      vq,
      te,
      id,
      iq
    };

    /**
     * @brief Contains modeling data for a Gensal generator model.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct GensalData : public ComponentData<real_type,
                                             index_type,
                                             GensalParameters,
                                             GensalBuses,
                                             GensalSignalInputs,
                                             GensalSignalOutputs,
                                             GensalMonitorableVariables>
    {
      GensalData() = default;

      using Parameters           = GensalParameters;
      using Buses                = GensalBuses;
      using SignalInputs         = GensalSignalInputs;
      using SignalOutputs        = GensalSignalOutputs;
      using MonitorableVariables = GensalMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
