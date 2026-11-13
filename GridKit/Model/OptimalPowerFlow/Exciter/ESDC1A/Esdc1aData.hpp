/**
 * @file Esdc1aData.hpp
 * @brief Modeling data for optimal power flow ESDC1A exciters.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      /// Parameters for an ESDC1A exciter, those its steady state reads
      enum class Esdc1aParameters : size_t
      {
        Ke,     ///< \f$K_E\f$ Exciter field resistance line slope margin [p.u.]
        Vrmax,  ///< \f$V_R^{\max}\f$ Maximum voltage-regulator output [p.u.]
        Vrmin,  ///< \f$V_R^{\min}\f$ Minimum voltage-regulator output [p.u.]
        E1,     ///< \f$E_1\f$ First saturation voltage point [p.u.]
        Se1,    ///< \f$S_E(E_1)\f$ Saturation coefficient at \f$E_1\f$ [p.u.]
        E2,     ///< \f$E_2\f$ Second saturation voltage point [p.u.]
        Se2,    ///< \f$S_E(E_2)\f$ Saturation coefficient at \f$E_2\f$ [p.u.]
        exclim, ///< \f$s_{\mathrm{lim}}\f$ Exciter field-voltage-state lower-limit flag
      };

      /// Signals of an ESDC1A exciter
      enum class Esdc1aSignals : size_t
      {
        efd, ///< Field voltage
      };

      /// Data of an `Esdc1a`
      template <typename real_type, typename index_type>
      using Esdc1aData = ComponentData<real_type, index_type, Esdc1aParameters, Utilities::EmptyEnum, Esdc1aSignals>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
