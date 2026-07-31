/**
 * @file Pf1Data.hpp
 * @brief Modeling data for the PF1 model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /// Parameter keys for the PF1 model.
      enum class Pf1Parameters
      {
        Tslew,     ///< \f$T_{\mathrm{slew}}\f$ Voltage adjuster travel time
        VREFmax,   ///< \f$V_{\mathrm{REFmax}}\f$ Voltage adjuster maximum output
        VREFmin,   ///< \f$V_{\mathrm{REFmin}}\f$ Voltage adjuster minimum output
        Ton,       ///< \f$T_{\mathrm{on}}\f$ Voltage adjuster pulse generator time on
        Toff,      ///< \f$T_{\mathrm{off}}\f$ Voltage adjuster pulse generator time off
        VADJF,     ///< \f$V_{\mathrm{ADJF}}\f$ Voltage adjuster bypass of pulse generator
        PFREFnorm, ///< \f$\mathrm{PF}_{\mathrm{REFnorm}}\f$ Power factor controller normalized reference setpoint
        VITmin,    ///< \f$V_{\mathrm{ITmin}}\f$ Power factor controller minimum terminal current limit
        VVTmin,    ///< \f$V_{\mathrm{VTmin}}\f$ Power factor controller minimum terminal voltage limit
        VVTmax,    ///< \f$V_{\mathrm{VTmax}}\f$ Power factor controller maximum terminal voltage limit
        VPFC_BW,   ///< \f$V_{\mathrm{PFC_BW}}\f$ Power factor controller deadband magnitude
        TPFC       ///< \f$T_{\mathrm{PFC}}\f$ Power factor controller delay time
      };

      /// Bus keys for the PF1 model; deferred until port integration.
      enum class Pf1Buses : size_t
      {
        SIZE
      };

      /// Signal input keys for the PF1 model; deferred until port integration.
      enum class Pf1SignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the PF1 model; deferred until port integration.
      enum class Pf1SignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the PF1 model; deferred until implementation.
      enum class Pf1MonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for PF1.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Pf1Data : public ComponentData<real_type,
                                            index_type,
                                            Pf1Parameters,
                                            Pf1Buses,
                                            Pf1SignalInputs,
                                            Pf1SignalOutputs,
                                            Pf1MonitorableVariables>
      {
        Pf1Data() = default;

        using Parameters           = Pf1Parameters;
        using Buses                = Pf1Buses;
        using SignalInputs         = Pf1SignalInputs;
        using SignalOutputs        = Pf1SignalOutputs;
        using MonitorableVariables = Pf1MonitorableVariables;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
