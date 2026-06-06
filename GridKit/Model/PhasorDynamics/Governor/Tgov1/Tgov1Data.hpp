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
       * @brief Temporary TGOV1 terminal keys.
       *
       * NOTE: The TGOV1 bus terminal is accepted only so existing flat JSON
       * cases continue to parse without case-file churn. TGOV1 does not use
       * this bus terminal, and it should be removed when the JSON port format
       * is updated.
       */
      enum class Tgov1Terminals : size_t
      {
        bus,
        SIZE,
      };

      /**
       * @brief TGOV1 input ports.
       */
      enum class Tgov1InputPorts : size_t
      {
        speed,
        SIZE,
      };

      /**
       * @brief TGOV1 output ports.
       */
      enum class Tgov1OutputPorts : size_t
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
      struct Tgov1Data : public ComponentData<real_type,
                                              index_type,
                                              Tgov1Parameters,
                                              Tgov1Terminals,
                                              Tgov1InputPorts,
                                              Tgov1OutputPorts,
                                              Tgov1MonitorableVariables>
      {
        Tgov1Data() = default;

        using Parameters           = Tgov1Parameters;
        using Terminals            = Tgov1Terminals;
        using InputPorts           = Tgov1InputPorts;
        using OutputPorts          = Tgov1OutputPorts;
        using MonitorableVariables = Tgov1MonitorableVariables;
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
