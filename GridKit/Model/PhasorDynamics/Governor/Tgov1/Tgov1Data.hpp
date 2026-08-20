/**
 * @file Tgov1Data.hpp
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for TGOV1
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Parameter keys for TGOV1 Governor model.
       *
       * These enum values serve as keys for the parameters map in ComponentData.
       */
      enum class Tgov1Parameters : size_t
      {
        Trate, ///< \f$T_\mathrm{rate}\f$ Turbine-rating power base
        R,     ///< \f$R\f$ Permanent droop
        T1,    ///< \f$T_1\f$ Steam-bowl time constant
        T2,    ///< \f$T_2\f$ Turbine numerator time constant
        T3,    ///< \f$T_3\f$ Reheater time constant
        Pvmax, ///< \f$P_v^\mathrm{max}\f$ Maximum valve position
        Pvmin, ///< \f$P_v^\mathrm{min}\f$ Minimum valve position
        Dt,    ///< \f$D_t\f$ Turbine damping coefficient
      };

      /**
       * @brief Placeholder enum for TGOV1 bus keys.
       */
      enum class Tgov1Buses : size_t
      {
      };

      /**
       * @brief TGOV1 signal inputs.
       */
      enum class Tgov1SignalInputs : size_t
      {
        speed, ///< \f$\omega\f$ Optional machine speed-deviation signal ID
        pref,  ///< \f$P_\mathrm{ref}\f$ Optional governor-reference signal ID
      };

      /**
       * @brief TGOV1 signal outputs.
       */
      enum class Tgov1SignalOutputs : size_t
      {
        pmech, ///< \f$P_m\f$ Required mechanical-power output signal ID
      };

      /**
       * @brief Placeholder enum for TGOV1 monitorable variables.
       */
      enum class Tgov1MonitorableVariables : size_t
      {
      };

      /**
       * @brief Modeling data for TGOV1 Governor using ComponentData base.
       *
       * @tparam real_type  Real number type (e.g., double)
       * @tparam index_type Index type (e.g., size_t)
       */
      template <typename real_type, typename index_type>
      using Tgov1Data =
          ComponentData<real_type,
                        index_type,
                        Tgov1Parameters,
                        Tgov1Buses,
                        Tgov1SignalInputs,
                        Tgov1SignalOutputs,
                        Tgov1MonitorableVariables>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
