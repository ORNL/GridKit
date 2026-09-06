/**
 * @file GenClassicalData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for a classical generator model.
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a classical generator model
    enum class GenClassicalParameters
    {
      p0,  ///< \f$P_0\f$ Initial active power
      q0,  ///< \f$Q_0\f$ Initial reactive power
      H,   ///< \f$H\f$ Rotor inertia
      D,   ///< \f$D\f$ Damping coefficient
      Ra,  ///< \f$R_a\f$ Winding resistance
      Xdp, ///< \f$X'_d\f$ Direct axis transient reactance
      mva, ///< \f$S_\mathrm{mach}\f$ MVA base of the classical generator model
    };

    /// Buses for a classical generator model
    enum class GenClassicalBuses : size_t
    {
      bus, ///< Unique ID of the connecting bus
      SIZE
    };

    /// Signal inputs for a classical generator model
    enum class GenClassicalSignalInputs : size_t
    {
      pmech, ///< \f$P_m\f$ Unique ID of the signal providing mechanical power
      efd,   ///< \f$E_{fd}\f$ Unique ID of the signal providing exciter field voltage
      SIZE
    };

    /// Signal outputs for a classical generator model
    enum class GenClassicalSignalOutputs : size_t
    {
      speed, ///< \f$\omega\f$ Unique ID of the signal receiving speed deviation
      SIZE
    };

    /// Variables able to be monitored for a classical generator model
    enum class GenClassicalMonitorableVariables
    {
      ir,    ///< \f$I_r\f$ Network-frame real terminal current
      ii,    ///< \f$I_i\f$ Network-frame imaginary terminal current
      p,     ///< \f$P\f$ Active power
      q,     ///< \f$Q\f$ Reactive power
      delta, ///< \f$\delta\f$ Rotor angle
      omega, ///< \f$\omega\f$ Speed deviation
      speed  ///< \f$1+\omega\f$ Per-unit machine speed
    };

    /**
     * @brief Contains modeling data for a classical generator model.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct GenClassicalData : public ComponentData<real_type,
                                                   index_type,
                                                   GenClassicalParameters,
                                                   GenClassicalBuses,
                                                   GenClassicalSignalInputs,
                                                   GenClassicalSignalOutputs,
                                                   GenClassicalMonitorableVariables>
    {
      GenClassicalData() = default;

      using Parameters           = GenClassicalParameters;
      using Buses                = GenClassicalBuses;
      using SignalInputs         = GenClassicalSignalInputs;
      using SignalOutputs        = GenClassicalSignalOutputs;
      using MonitorableVariables = GenClassicalMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
