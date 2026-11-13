/**
 * @file SexsPtiData.hpp
 * @brief Modeling data for optimal power flow SEXS-PTI exciters.
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
      /// Parameters for a SEXS-PTI exciter, those its steady state reads
      enum class SexsPtiParameters : size_t
      {
        Efdmin, ///< \f$E_\mathrm{fd}^{\min}\f$ Minimum excitation output [p.u.]
        Efdmax, ///< \f$E_\mathrm{fd}^{\max}\f$ Maximum excitation output [p.u.]
      };

      /// Signals of a SEXS-PTI exciter
      enum class SexsPtiSignals : size_t
      {
        efd, ///< Field voltage
      };

      /// Data of a `SexsPti`
      template <typename real_type, typename index_type>
      using SexsPtiData = ComponentData<real_type, index_type, SexsPtiParameters, Utilities::EmptyEnum, SexsPtiSignals>;
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
