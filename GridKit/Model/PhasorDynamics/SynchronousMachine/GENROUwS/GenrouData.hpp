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
      p0,       ///< Initial active power
      q0,       ///< Initial reactive power
      H,        ///< Rotor inertia
      D,        ///< Damping coefficient
      Ra,       ///< Winding resistance
      Tdop,     ///< Open circuit direct axis transient time
      Tdopp,    ///< Open circuit direct axis sub-transient time
      Tqop,     ///< Open circuit quadrature axis transient
      Tqopp,    ///< Open circuit quadrature axis sub-transient time
      Xd,       ///< Direct axis synchronous reactance
      Xdp,      ///< Direct axis transient reactance
      Xdpp,     ///< Direct axis sub-transient reactance
      Xq,       ///< Quadrature axis synchronous reactance
      Xqp,      ///< Quadrature axis transient reactance
      Xqpp,     ///< Quadrature axis sub-transient reactance
      Xl,       ///< Stator leakage reactance
      S10,      ///< Saturation factor at 1.0 pu flux
      S12,      ///< Saturation factor at 1.2 pu flux
      mva_base, ///< MVA base of the genrou model
    };

    /// Ports for a Genrou generator model
    enum class GenrouPorts
    {
      bus,   ///< Unique ID of the connecting bus
      pmech, ///< Unique ID of the bus providing the exciter signal
      speed, ///< Unique ID of the bus providing the governor signal
      efd,   ///< Unique ID of the bus providing exciter field signal
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
