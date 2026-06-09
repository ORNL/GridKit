/**
 * @file HygovData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the HYGOV governor model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /// Parameter keys for the HYGOV governor model.
      enum class HygovParameters
      {
        Trate, ///< Turbine-rating power base
        mva,   ///< MVA base of the pmech signal
        Rperm, ///< Permanent droop
        Rtemp, ///< Temporary droop
        Tr,    ///< Temporary-droop reset time constant
        Tf,    ///< Governor error filter time constant
        Tg,    ///< Gate servo time constant
        Velm,  ///< Maximum desired-gate velocity magnitude
        Gmax,  ///< Maximum desired-gate position
        Gmin,  ///< Minimum desired-gate position
        Tw,    ///< Water inertia time constant
        At,    ///< Turbine gain
        Dturb, ///< Turbine damping coefficient
        Qnl,   ///< No-load flow at nominal head
        Tn,    ///< Speed lead-lag numerator time constant
        Tnp,   ///< Speed lead-lag denominator time constant
        db1,   ///< Type 1 speed deadband threshold
        db2,   ///< Unsupported mechanical backlash deadband
        Hdam,  ///< Head available at dam
        Gv0,   ///< Gate point 0
        Gv1,   ///< Gate point 1
        Gv2,   ///< Gate point 2
        Gv3,   ///< Gate point 3
        Gv4,   ///< Gate point 4
        Gv5,   ///< Gate point 5
        Pgv0,  ///< Power point 0
        Pgv1,  ///< Power point 1
        Pgv2,  ///< Power point 2
        Pgv3,  ///< Power point 3
        Pgv4,  ///< Power point 4
        Pgv5   ///< Power point 5
      };

      /// Ports for the HYGOV governor model.
      enum class HygovPorts
      {
        bus,   ///< Optional terminal bus ID for case-format compatibility
        speed, ///< Machine speed-deviation signal ID
        pmech, ///< Mechanical-power output signal ID
        pref,  ///< Optional active-power/load reference signal ID
        paux   ///< Optional auxiliary power input signal ID
      };

      /// Variables available through the monitor interface.
      enum class HygovMonitorableVariables
      {
        pmech,       ///< Mechanical power output
        leadlag,     ///< Speed lead-lag denominator state
        filter,      ///< Governor error filter output
        desiredgate, ///< Desired-gate position
        gate,        ///< Gate position
        flow,        ///< Turbine flow
        head,        ///< Turbine head
        pgv          ///< Nonlinear gate-to-power curve output
      };

      template <typename real_type, typename index_type>
      struct HygovData : public ComponentData<real_type,
                                              index_type,
                                              HygovParameters,
                                              HygovPorts,
                                              HygovMonitorableVariables>
      {
        HygovData() = default;

        using Parameters           = HygovParameters;
        using Ports                = HygovPorts;
        using MonitorableVariables = HygovMonitorableVariables;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
