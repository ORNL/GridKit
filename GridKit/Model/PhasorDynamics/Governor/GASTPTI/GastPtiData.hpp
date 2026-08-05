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
      /// Parameter keys for GASTPTI; all are optional with documented defaults.
      enum class GastPtiParameters
      {
        R,     ///< \f$R\f$ Speed-deviation droop per component-power deviation [p.u.]
        T1,    ///< \f$T_1\f$ Fuel-valve time constant [sec]
        T2,    ///< \f$T_2\f$ Fuel-flow time constant [sec]
        T3,    ///< \f$T_3\f$ Exhaust-temperature time constant [sec]
        At,    ///< \f$A_T\f$ Ambient-temperature load limit on component base [p.u.]
        Kt,    ///< \f$K_T\f$ Exhaust-temperature feedback gain [p.u.]
        Vmax,  ///< \f$V^{\max}\f$ Configured upper valve limit on component base [p.u.]
        Vmin,  ///< \f$V^{\min}\f$ Configured lower valve limit on component base [p.u.]
        Dturb, ///< \f$D^\mathrm{turb}\f$ Component-base power per speed deviation [p.u.]
        Trate  ///< \f$T^\mathrm{rate}\f$ MW rating used as the same-valued MVA component base [MW]
      };

      /// Buses for the GASTPTI governor model.
      enum class GastPtiBuses : size_t
      {
        SIZE ///< Number of GASTPTI bus ports; GASTPTI has no bus attachment
      };

      /// Signal inputs for the GASTPTI governor model.
      enum class GastPtiSignalInputs : size_t
      {
        speed, ///< \f$\omega\f$ Optional Known machine speed-deviation input [p.u.]
        pref,  ///< \f$P^\mathrm{ref}\f$ Optional Unknown load-reference input on system base [p.u.]
        SIZE   ///< Number of GASTPTI signal-input ports
      };

      /// Signal outputs for the GASTPTI governor model.
      enum class GastPtiSignalOutputs : size_t
      {
        pmech, ///< \f$P_{\text{m}}\f$ Required Known mechanical-power output on system base [p.u.]
        SIZE   ///< Number of GASTPTI signal-output ports
      };

      /// Variables available through the monitor interface.
      enum class GastPtiMonitorableVariables
      {
        pmech,  ///< \f$P_{\text{m}}\f$ Mechanical-power output on system base [p.u.]
        xvalve, ///< \f$x_V\f$ Fuel-valve state on component base [p.u.]
        xflow,  ///< \f$x_F\f$ Fuel-flow state on component base [p.u.]
        xtemp,  ///< \f$x_T\f$ Exhaust-temperature feedback state on component base [p.u.]
        vload,  ///< \f$V_D\f$ Speed/load fuel demand on component base [p.u.]
        vtemp   ///< \f$V_T\f$ Temperature-limit fuel demand on component base [p.u.]
      };

      /**
       * @brief Model data for GASTPTI parameters, signal ports, and monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see GastPti
       */
      template <typename real_type, typename index_type>
      struct GastPtiData : public ComponentData<real_type,
                                                index_type,
                                                GastPtiParameters,
                                                GastPtiBuses,
                                                GastPtiSignalInputs,
                                                GastPtiSignalOutputs,
                                                GastPtiMonitorableVariables>
      {
        GastPtiData() = default;

        using Parameters           = GastPtiParameters;
        using Buses                = GastPtiBuses;
        using SignalInputs         = GastPtiSignalInputs;
        using SignalOutputs        = GastPtiSignalOutputs;
        using MonitorableVariables = GastPtiMonitorableVariables;
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
