/**
 * @file RepcaData.hpp
 * @brief Modeling data for optimal power flow REPCA controllers.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      /// Parameters for a REPCA controller, those its steady state reads
      enum class RepcaParameters : size_t
      {
        mva,      ///< \f$S^\mathrm{base}\f$ REPCA component power base [MVA]
        Freqflag, ///< \f$s_\mathrm{freq}\f$ Active-power output selector
        Qmax,     ///< \f$Q^{\max}\f$ Maximum reactive-power command [p.u.]
        Qmin,     ///< \f$Q^{\min}\f$ Minimum reactive-power command [p.u.]
        Pmax,     ///< \f$P^{\max}\f$ Maximum active-power command [p.u.]
        Pmin,     ///< \f$P^{\min}\f$ Minimum active-power command [p.u.]
      };

      /// Signals of a REPCA controller
      enum class RepcaSignals : size_t
      {
        p, ///< Branch active power
        q, ///< Branch reactive power
      };

      /// Data of a `Repca`
      template <typename real_type, typename index_type>
      using RepcaData = ComponentData<real_type, index_type, RepcaParameters, Utilities::EmptyEnum, RepcaSignals>;
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
