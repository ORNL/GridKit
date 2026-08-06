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
       * @brief Placeholder enum for TGOV1 bus keys.
       */
      enum class Tgov1Buses : size_t
      {
        SIZE,
      };

      /**
       * @brief TGOV1 signal inputs.
       */
      enum class Tgov1SignalInputs : size_t
      {
        speed, ///< Optional machine speed-deviation signal ID
        pref,  ///< Optional governor-reference signal ID
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
      struct Tgov1Data : public ComponentData<real_type,
                                              index_type,
                                              Tgov1Parameters,
                                              Tgov1Buses,
                                              Tgov1SignalInputs,
                                              Tgov1SignalOutputs,
                                              Tgov1MonitorableVariables>
      {
        Tgov1Data() = default;

        using Parameters           = Tgov1Parameters;
        using Buses                = Tgov1Buses;
        using SignalInputs         = Tgov1SignalInputs;
        using SignalOutputs        = Tgov1SignalOutputs;
        using MonitorableVariables = Tgov1MonitorableVariables;
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
