/**
 * @file Tgov1Data.hpp
 * @author Wiktoria Zielinska (zielinskawa@ORNL.gov)
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for TGOV1
 */

#pragma once

#include <Model/PhasorDynamics/ComponentData.hpp>

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
        R,     ///< Droop Constant
        T1,    ///< Valve Time Delay
        T2,    ///< Turbine Numerator Time Constant
        T3,    ///< Turbine Delay
        Pvmax, ///< Max Valve Power
        Pvmin, ///< Min Valve Power
        Dt     ///< Damping Coefficient
      };

      /**
       * @brief Placeholder enum for TGOV1 ports.
       */
      enum class Tgov1Ports
      {
        signal_in,
        signal_out,
      };

      /**
       * @brief Placeholder enum for TGOV1 monitorable variables.
       */
      enum class Tgov1MonitorableVariables
      {
      };

      /**
       * @brief Modeling data for TGOV1 Governor using ComponentData base.
       *
       * @tparam RealT Real number type (e.g., double)
       * @tparam IdxT  Index type (e.g., size_t)
       */
      template <typename RealT, typename IdxT>
      struct Tgov1Data : public ComponentData<RealT, IdxT, Tgov1Parameters, Tgov1Ports, Tgov1MonitorableVariables>
      {
        Tgov1Data() = default;

        using Parameters           = Tgov1Parameters;
        using Ports                = Tgov1Ports;
        using MonitorableVariables = Tgov1MonitorableVariables;
      };

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
