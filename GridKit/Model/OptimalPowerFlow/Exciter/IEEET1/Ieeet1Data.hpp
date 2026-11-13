/**
 * @file Ieeet1Data.hpp
 * @brief Modeling data for optimal power flow IEEET1 exciters.
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
      /// Parameters for an IEEET1 exciter, those its steady state reads
      enum class Ieeet1Parameters : size_t
      {
        Ke,    ///< \f$K_E\f$ Exciter field resistance line slope margin [p.u.]
        Vrmin, ///< \f$V_R^{\min}\f$ Lower limit to voltage regulation [p.u.]
        Vrmax, ///< \f$V_R^{\max}\f$ Upper limit to voltage regulation [p.u.]
        E1,    ///< \f$E_1\f$ Saturation parameter [p.u.]
        E2,    ///< \f$E_2\f$ Saturation parameter [p.u.]
        Se1,   ///< \f$S_1\f$ Saturation parameter [p.u.]
        Se2,   ///< \f$S_2\f$ Saturation parameter [p.u.]
      };

      /// Signals of an IEEET1 exciter
      enum class Ieeet1Signals : size_t
      {
        efd, ///< Field voltage
      };

      /// Data of an `Ieeet1`
      template <typename real_type, typename index_type>
      using Ieeet1Data = ComponentData<real_type, index_type, Ieeet1Parameters, Utilities::EmptyEnum, Ieeet1Signals>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
