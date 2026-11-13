/**
 * @file HygovData.hpp
 * @brief Modeling data for optimal power flow HYGOV governors.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      /// Parameters for a HYGOV governor, those its steady state reads
      enum class HygovParameters : size_t
      {
        Trate, ///< \f$T^\mathrm{rate}\f$ Turbine-rating power base [MW]
        Gmax,  ///< \f$G^{\max}\f$ Configured upper gate response limit [p.u.]
        Gmin,  ///< \f$G^{\min}\f$ Configured lower gate response limit [p.u.]
        At,    ///< \f$A_\mathrm{t}\f$ Turbine gain [p.u.]
        Qnl,   ///< \f$q_{\mathrm{NL}}\f$ No-load flow at nominal head [p.u.]
        Hdam,  ///< \f$H_{\mathrm{dam}}\f$ Configured dam head [p.u.]
      };

      /// Signals of a HYGOV governor
      enum class HygovSignals : size_t
      {
        pmech, ///< Mechanical power
      };

      /// Data of a `Hygov`
      template <typename real_type, typename index_type>
      using HygovData = ComponentData<real_type, index_type, HygovParameters, Utilities::EmptyEnum, HygovSignals>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
