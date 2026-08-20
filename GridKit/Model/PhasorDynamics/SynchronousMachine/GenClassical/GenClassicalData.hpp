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
    enum class GenClassicalParameters : size_t
    {
      p0,  ///< Initial active power
      q0,  ///< Initial reactive power
      H,   ///< Rotor inertia
      D,   ///< Damping coefficient
      Ra,  ///< Winding resistance
      Xdp, ///< Direct axis transient reactance
      mva, ///< MVA Base of the generator
    };

    /// Buses supported for a classical generator model
    enum class GenClassicalBuses : size_t
    {
      bus, ///< Unique ID of the connecting bus
    };

    /// Signal inputs supported for a classical generator model
    ///
    /// @warning GenClassical signal support is incomplete. These legacy signal
    /// names are not wired by SystemModel today; the intended refactor is to
    /// align this model with Genrou/Gensal by supporting `pmech`, `speed`, and
    /// `efd` signals.
    enum class GenClassicalSignalInputs : size_t
    {
      exciter_signal,  ///< Unique ID of the bus providing the exciter signal
      governor_signal, ///< Unique ID of the bus providing the governor signal
    };

    /// Signal outputs supported for a classical generator model
    enum class GenClassicalSignalOutputs : size_t
    {
    };

    /// Variables able to be monitored for a classical generator model
    enum class GenClassicalMonitorableVariables : size_t
    {
      ir,
      ii,
      p,
      q,
      delta,
      omega,
      speed,
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
    using GenClassicalData =
        ComponentData<real_type,
                      index_type,
                      GenClassicalParameters,
                      GenClassicalBuses,
                      GenClassicalSignalInputs,
                      GenClassicalSignalOutputs,
                      GenClassicalMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
