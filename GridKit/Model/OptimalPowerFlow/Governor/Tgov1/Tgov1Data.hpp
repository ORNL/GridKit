/**
 * @file Tgov1Data.hpp
 * @brief Modeling data for optimal power flow TGOV1 governors.
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
      /// Parameters for a TGOV1 governor, those its steady state reads
      enum class Tgov1Parameters : size_t
      {
        Trate, ///< \f$T_{\mathrm{rate}}\f$ Governor component power base [MVA]
        Pvmax, ///< \f$P_\mathrm{v}^{\max}\f$ Maximum valve position [p.u.]
        Pvmin, ///< \f$P_\mathrm{v}^{\min}\f$ Minimum valve position [p.u.]
      };

      /// Signals of a TGOV1 governor
      enum class Tgov1Signals : size_t
      {
        pmech, ///< Mechanical power
      };

      /// Data of a `Tgov1`
      template <typename real_type, typename index_type>
      using Tgov1Data = ComponentData<real_type, index_type, Tgov1Parameters, Utilities::EmptyEnum, Tgov1Signals>;
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
