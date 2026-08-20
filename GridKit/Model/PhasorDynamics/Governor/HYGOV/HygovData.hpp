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
      enum class HygovParameters : size_t
      {
        Trate, ///< \f$T^\mathrm{rate}\f$ Turbine-rating power base
        Rperm, ///< \f$R_{\mathrm{perm}}\f$ Permanent droop
        Rtemp, ///< \f$R_{\mathrm{temp}}\f$ Temporary droop
        Tr,    ///< \f$T_r\f$ Temporary-droop reset time constant
        Tf,    ///< \f$T_f\f$ Governor error filter time constant
        Tg,    ///< \f$T_g\f$ Gate servo time constant
        Velm,  ///< \f$V_{\mathrm{elm}}\f$ Maximum desired-gate velocity magnitude
        Gmax,  ///< \f$G^{\max}\f$ Configured upper gate response limit
        Gmin,  ///< \f$G^{\min}\f$ Configured lower gate response limit
        Tw,    ///< \f$T_w\f$ Water inertia time constant
        At,    ///< \f$A_t\f$ Turbine gain
        Dturb, ///< \f$D_{\mathrm{turb}}\f$ Turbine damping coefficient
        Qnl,   ///< \f$q_{\mathrm{NL}}\f$ No-load flow at nominal head
        Tn,    ///< \f$T_n\f$ Speed lead-lag numerator time constant
        Tnp,   ///< \f$T_{\mathrm{np}}\f$ Speed lead-lag denominator time constant
        db1,   ///< \f$D_{\omega}\f$ Type 1 speed deadband threshold
        db2,   ///< \f$D_2\f$ Unsupported mechanical backlash. Nonzero values warn and are ignored
        Hdam,  ///< \f$H_{\mathrm{dam}}\f$ Head available at dam
        Gv0,   ///< \f$G_V^{(0)}\f$ Gate point 0
        Gv1,   ///< \f$G_V^{(1)}\f$ Gate point 1
        Gv2,   ///< \f$G_V^{(2)}\f$ Gate point 2
        Gv3,   ///< \f$G_V^{(3)}\f$ Gate point 3
        Gv4,   ///< \f$G_V^{(4)}\f$ Gate point 4
        Gv5,   ///< \f$G_V^{(5)}\f$ Gate point 5
        Pgv0,  ///< \f$P_{\mathrm{GV}}^{(0)}\f$ Power point 0
        Pgv1,  ///< \f$P_{\mathrm{GV}}^{(1)}\f$ Power point 1
        Pgv2,  ///< \f$P_{\mathrm{GV}}^{(2)}\f$ Power point 2
        Pgv3,  ///< \f$P_{\mathrm{GV}}^{(3)}\f$ Power point 3
        Pgv4,  ///< \f$P_{\mathrm{GV}}^{(4)}\f$ Power point 4
        Pgv5,  ///< \f$P_{\mathrm{GV}}^{(5)}\f$ Power point 5
      };

      /// Buses for the HYGOV governor model.
      enum class HygovBuses : size_t
      {
      };

      /// Signal inputs for the HYGOV governor model.
      enum class HygovSignalInputs : size_t
      {
        speed, ///< Optional machine speed-deviation signal ID
        pref,  ///< Optional active-power/load reference signal ID
        paux,  ///< Optional auxiliary power input signal ID
      };

      /// Signal outputs for the HYGOV governor model.
      enum class HygovSignalOutputs : size_t
      {
        pmech, ///< Required mechanical-power output signal ID
      };

      /// Variables available through the monitor interface.
      enum class HygovMonitorableVariables : size_t
      {
        pmech,       ///< Mechanical power output on system base
        filter,      ///< Governor error filter output on component base
        desiredgate, ///< Desired-gate position on component base
        gate,        ///< Gate position on component base
        flow,        ///< Turbine flow on component base
        head,        ///< Turbine head on component base
      };

      /**
       * @brief Model data for HYGOV: parameters, optional input signals, the
       *        required mechanical-power output, and monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see Hygov
       */
      template <typename real_type, typename index_type>
      struct HygovData : public ComponentData<real_type,
                                              index_type,
                                              HygovParameters,
                                              HygovBuses,
                                              HygovSignalInputs,
                                              HygovSignalOutputs,
                                              HygovMonitorableVariables>
      {
        HygovData() = default;

        using Parameters           = HygovParameters;
        using Buses                = HygovBuses;
        using SignalInputs         = HygovSignalInputs;
        using SignalOutputs        = HygovSignalOutputs;
        using MonitorableVariables = HygovMonitorableVariables;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
