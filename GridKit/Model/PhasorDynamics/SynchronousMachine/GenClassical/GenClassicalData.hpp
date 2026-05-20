/**
 * @file GenClassicalData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
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
      p0,  ///< Initial active power
      q0,  ///< Initial reactive power
      H,   ///< Rotor inertia
      D,   ///< Damping coefficient
      Ra,  ///< Winding resistance
      Xdp, ///< Direct axis transient reactance
      mva  ///< MVA Base of the generator
    };

    /// Ports supported for a classical generator model
    ///
    /// @warning GenClassical signal support is incomplete. These legacy port
    /// names are not wired by SystemModel today; the intended refactor is to
    /// align this model with Genrou/Gensal by supporting `pmech`, `speed`, and
    /// `efd` ports through ComponentSignals.
    enum class GenClassicalPorts
    {
      bus,             ///< Unique ID of the connecting bus
      exciter_signal,  ///< Unique ID of the bus providing the exciter signal
      governor_signal, ///< Unique ID of the bus providing the governor signal
    };

    /// Variables able to be monitored for a classical generator model
    enum class GenClassicalMonitorableVariables
    {
      ir,
      ii,
      p,
      q,
      delta,
      omega
    };

    /**
     * @brief Contains modeling data for a GenClassical generator model.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename real_type, typename index_type>
    struct GenClassicalData : public ComponentData<real_type,
                                                   index_type,
                                                   GenClassicalParameters,
                                                   GenClassicalPorts,
                                                   GenClassicalMonitorableVariables>
    {
      GenClassicalData() = default;

      using Parameters           = GenClassicalParameters;
      using Ports                = GenClassicalPorts;
      using MonitorableVariables = GenClassicalMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
