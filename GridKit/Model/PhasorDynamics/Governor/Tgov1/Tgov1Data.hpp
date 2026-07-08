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
      enum class Tgov1Parameters
      {
        Trate, ///< Turbine-rating power base
        R,     ///< Droop Constant
        T1,    ///< Valve Time Delay
        T2,    ///< Turbine Numerator Time Constant
        T3,    ///< Turbine Delay
        Pvmax, ///< Max Valve Power
        Pvmin, ///< Min Valve Power
        Dt     ///< Damping Coefficient
      };

      /**
       * @brief Temporary TGOV1 bus keys.
       *
       * NOTE: The TGOV1 bus key is accepted only so existing flat JSON
       * cases continue to parse without case-file churn. TGOV1 does not use
       * this bus key, and it should be removed when the JSON port format
       * is updated.
       */
      enum class Tgov1Buses : size_t
      {
        bus,
        SIZE,
      };

      /**
       * @brief TGOV1 signal inputs.
       */
      enum class Tgov1SignalInputs : size_t
      {
        speed,
        SIZE,
      };

      /**
       * @brief TGOV1 signal outputs.
       */
      enum class Tgov1SignalOutputs : size_t
      {
        pmech,
        SIZE,
      };

      /**
       * @brief Placeholder enum for TGOV1 monitorable variables.
       */
      enum class Tgov1MonitorableVariables
      {
        NONE,
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
