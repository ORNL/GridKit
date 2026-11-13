/**
 * @file GensalData.hpp
 * @brief Modeling data for optimal power flow salient-pole machines.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a salient-pole machine, those its steady state reads
    enum class GensalParameters : size_t
    {
      Ra,  ///< \f$R_\mathrm{a}\f$ Winding resistance [p.u.]
      Xd,  ///< \f$X_d\f$ Direct axis synchronous reactance [p.u.]
      Xdp, ///< \f$X'_d\f$ Direct axis transient reactance [p.u.]
      Xq,  ///< \f$X_q\f$ Quadrature axis synchronous reactance [p.u.]
      S10, ///< \f$S_{10}\f$ Saturation factor at 1.0 pu flux [p.u.]
      S12, ///< \f$S_{12}\f$ Saturation factor at 1.2 pu flux [p.u.]
      mva, ///< \f$S^\mathrm{base}\f$ Machine power base [MVA]
    };

    /// Ports of a salient-pole machine
    enum class GensalPorts : size_t
    {
      bus, ///< Port bus
    };

    /// Signals of a salient-pole machine
    enum class GensalSignals : size_t
    {
      pmech, ///< Mechanical power
      efd,   ///< Field voltage
    };

    /// Data of a `Gensal`
    template <typename real_type, typename index_type>
    using GensalData = ComponentData<real_type, index_type, GensalParameters, GensalPorts, GensalSignals>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
