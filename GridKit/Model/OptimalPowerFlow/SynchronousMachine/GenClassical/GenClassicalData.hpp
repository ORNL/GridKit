/**
 * @file GenClassicalData.hpp
 * @brief Modeling data for optimal power flow classical machines.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a classical machine, those its steady state reads
    enum class GenClassicalParameters : size_t
    {
      Ra,  ///< \f$R_\mathrm{a}\f$ Armature resistance [p.u.]
      Xdp, ///< \f$X'_d\f$ Direct-axis transient reactance [p.u.]
      mva, ///< \f$S^\mathrm{base}\f$ Machine power base [MVA]
    };

    /// Ports of a classical machine
    enum class GenClassicalPorts : size_t
    {
      bus, ///< Port bus
    };

    /// Signals of a classical machine
    enum class GenClassicalSignals : size_t
    {
      pmech, ///< Mechanical power
      efd,   ///< Field voltage
    };

    /// Data of a `GenClassical`
    template <typename real_type, typename index_type>
    using GenClassicalData = ComponentData<real_type, index_type, GenClassicalParameters, GenClassicalPorts, GenClassicalSignals>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
