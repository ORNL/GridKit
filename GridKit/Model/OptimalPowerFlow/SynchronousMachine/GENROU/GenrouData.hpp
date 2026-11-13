/**
 * @file GenrouData.hpp
 * @brief Modeling data for optimal power flow round-rotor machines.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a round-rotor machine, those its steady state reads
    enum class GenrouParameters : size_t
    {
      Ra,   ///< \f$R_\mathrm{a}\f$ Winding resistance [p.u.]
      Xd,   ///< \f$X_d\f$ Direct axis synchronous reactance [p.u.]
      Xq,   ///< \f$X_q\f$ Quadrature axis synchronous reactance [p.u.]
      Xdpp, ///< \f$X''_d\f$ Direct axis sub-transient reactance [p.u.]
      Xl,   ///< \f$X_{\ell}\f$ Stator leakage reactance [p.u.]
      S10,  ///< \f$S_{10}\f$ Saturation factor at 1.0 pu flux [p.u.]
      S12,  ///< \f$S_{12}\f$ Saturation factor at 1.2 pu flux [p.u.]
      mva,  ///< \f$S^\mathrm{base}\f$ Machine power base [MVA]
    };

    /// Ports of a round-rotor machine
    enum class GenrouPorts : size_t
    {
      bus, ///< Port bus
    };

    /// Signals of a round-rotor machine
    enum class GenrouSignals : size_t
    {
      pmech, ///< Mechanical power
      efd,   ///< Field voltage
    };

    /// Data of a `Genrou`
    template <typename real_type, typename index_type>
    using GenrouData = ComponentData<real_type, index_type, GenrouParameters, GenrouPorts, GenrouSignals>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
