/**
 * @file ReecbData.hpp
 * @brief Modeling data for optimal power flow REECB controllers.
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
      /// Parameters for a REECB controller, those its steady state reads
      enum class ReecbParameters : size_t
      {
        mva,  ///< \f$S^\mathrm{base}\f$ REECB component power base [MVA]
        Pmax, ///< \f$P^{\max}\f$ Maximum active-power order [p.u.]
        Pmin, ///< \f$P^{\min}\f$ Minimum active-power order [p.u.]
        Qmax, ///< \f$Q^{\max}\f$ Maximum reactive-power control output [p.u.]
        Qmin, ///< \f$Q^{\min}\f$ Minimum reactive-power control output [p.u.]
        Imax, ///< \f$I^{\max}\f$ Maximum converter current [p.u.]
      };

      /// Signals of a REECB controller
      enum class ReecbSignals : size_t
      {
        ipcmd, ///< Active-current command
        iqcmd, ///< Reactive-current command
        pe,    ///< Active-power feedback
        qgen,  ///< Reactive-power feedback
      };

      /// Data of a `Reecb`
      template <typename real_type, typename index_type>
      using ReecbData = ComponentData<real_type, index_type, ReecbParameters, Utilities::EmptyEnum, ReecbSignals>;
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
