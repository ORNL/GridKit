/**
 * @file GenrouData.hpp
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
    /// Initial parameters for a Genrou generator model
    enum class GenrouParameters
    {
      p0,    ///< Initial active power
      q0,    ///< Initial reactive power
      H,     ///< Rotor inertia
      D,     ///< Damping coefficient
      Ra,    ///< Winding resistance
      Tdop,  ///< Open circuit direct axis transient time
      Tdopp, ///< Open circuit direct axis sub-transient time
      Tqop,  ///< Open circuit quadrature axis transient
      Tqopp, ///< Open circuit quadrature axis sub-transient time
      Xd,    ///< Direct axis synchronous reactance
      Xdp,   ///< Direct axis transient reactance
      Xdpp,  ///< Direct axis sub-transient reactance
      Xq,    ///< Quadrature axis synchronous reactance
      Xqp,   ///< Quadrature axis transient reactance
      Xqpp,  ///< Quadrature axis sub-transient reactance
      Xl,    ///< Stator leakage reactance
      S10,   ///< Saturation factor at 1.0 pu flux
      S12,   ///< Saturation factor at 1.2 pu flux
      mva,   ///< MVA base of the genrou model
    };

    /// Buses for a Genrou generator model
    enum class GenrouBuses : size_t
    {
      bus, ///< Unique ID of the connecting bus
      SIZE
    };

    /// Signal inputs for a Genrou generator model
    ///
    /// @todo Genrou (and likely other components) would need to name multiple
    /// signal inlets and outlets. For now we have only speed out and
    /// mechanical power in.
    enum class GenrouSignalInputs : size_t
    {
      pmech, ///< Unique ID of the signal providing mechanical power
      efd,   ///< Unique ID of the signal providing exciter field signal
      SIZE
    };

    /// Signal outputs for a Genrou generator model
    enum class GenrouSignalOutputs : size_t
    {
      speed, ///< Unique ID of the signal receiving speed deviation
      SIZE
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
      speed
    };

    /**
     * @brief Contains modeling data for a Genrou generator model.
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    using GenrouData =
        ComponentData<real_type,
                      index_type,
                      GenrouParameters,
                      GenrouBuses,
                      GenrouSignalInputs,
                      GenrouSignalOutputs,
                      GenrouMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
