/**
 * @file GastPtiData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the GASTPTI governor model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /// Parameter keys for the GASTPTI governor model.
      enum class GastPtiParameters
      {
        R,     ///< Permanent droop
        T1,    ///< Fuel-valve time constant
        T2,    ///< Fuel-flow time constant
        T3,    ///< Exhaust-temperature time constant
        At,    ///< Ambient-temperature load limit
        Kt,    ///< Exhaust-temperature feedback gain
        Vmax,  ///< Maximum fuel-valve/turbine-power limit
        Vmin,  ///< Minimum fuel-valve/turbine-power limit
        Dturb, ///< Turbine damping coefficient
        Trate  ///< Turbine-rating power base
      };

      /// Ports for the GASTPTI governor model.
      enum class GastPtiPorts
      {
        bus,   ///< Optional terminal bus ID for case-format compatibility
        speed, ///< Machine speed-deviation signal ID
        pmech, ///< Mechanical-power output signal ID
        pref   ///< Optional active-power/load reference signal ID
      };

      /// Variables available through the monitor interface.
      enum class GastPtiMonitorableVariables
      {
        pmech,       ///< Mechanical power output
        fuelvalve,   ///< Fuel-valve state
        fuelflow,    ///< Fuel-flow state
        exhausttemp, ///< Exhaust-temperature feedback state
        vload,       ///< Speed/load fuel demand
        vtemp,       ///< Temperature-limit fuel demand
        vlv          ///< Low-value gate output
      };

      template <typename real_type, typename index_type>
      struct GastPtiData : public ComponentData<real_type,
                                                index_type,
                                                GastPtiParameters,
                                                GastPtiPorts,
                                                GastPtiMonitorableVariables>
      {
        GastPtiData() = default;

        using Parameters           = GastPtiParameters;
        using Ports                = GastPtiPorts;
        using MonitorableVariables = GastPtiMonitorableVariables;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
