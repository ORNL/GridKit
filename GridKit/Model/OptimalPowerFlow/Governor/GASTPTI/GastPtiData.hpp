/**
 * @file GastPtiData.hpp
 * @brief Modeling data for optimal power flow GASTPTI governors.
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
      /// Parameters for a GASTPTI governor, those its steady state reads
      enum class GastPtiParameters : size_t
      {
        Trate, ///< \f$T^\mathrm{rate}\f$ Turbine rating [MW]
        Vmax,  ///< \f$V^{\max}\f$ Upper valve response limit [p.u.]
        Vmin,  ///< \f$V^{\min}\f$ Lower valve response limit [p.u.]
        At,    ///< \f$A_T\f$ Ambient-temperature load limit [p.u.]
      };

      /// Signals of a GASTPTI governor
      enum class GastPtiSignals : size_t
      {
        pmech, ///< Mechanical power
      };

      /// Data of a `GastPti`
      template <typename real_type, typename index_type>
      using GastPtiData = ComponentData<real_type, index_type, GastPtiParameters, Utilities::EmptyEnum, GastPtiSignals>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
