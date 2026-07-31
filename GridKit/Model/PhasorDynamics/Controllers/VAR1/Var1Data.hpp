/**
 * @file Var1Data.hpp
 * @brief Modeling data for the VAR1 model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Controller
    {
      /// Parameter keys for the VAR1 model.
      enum class Var1Parameters
      {
        Tslew,    ///< \f$T_{\mathrm{slew}}\f$ Voltage adjuster travel time
        VREFmax,  ///< \f$V_{\mathrm{REFmax}}\f$ Voltage adjuster maximum output
        VREFmin,  ///< \f$V_{\mathrm{REFmin}}\f$ Voltage adjuster minimum output
        Ton,      ///< \f$T_{\mathrm{on}}\f$ Voltage adjuster pulse generator time on
        Toff,     ///< \f$T_{\mathrm{off}}\f$ Voltage adjuster pulse generator time off
        VADJF,    ///< \f$V_{\mathrm{ADJF}}\f$ Voltage adjuster bypass of pulse generator
        QREF,     ///< \f$Q_{\mathrm{REF}}\f$ Var controller reference setpoint
        VITmin,   ///< \f$V_{\mathrm{ITmin}}\f$ Var controller minimum terminal current limit
        VVTmin,   ///< \f$V_{\mathrm{VTmin}}\f$ Var controller minimum terminal voltage limit
        VVTmax,   ///< \f$V_{\mathrm{VTmax}}\f$ Var controller maximum terminal voltage limit
        VVARC_BW, ///< \f$V_{\mathrm{VARC_BW}}\f$ Var controller deadband magnitude
        TVARC     ///< \f$T_{\mathrm{VARC}}\f$ Var controller delay time
      };

      /// Bus keys for the VAR1 model; deferred until port integration.
      enum class Var1Buses : size_t
      {
        SIZE
      };

      /// Signal input keys for the VAR1 model; deferred until port integration.
      enum class Var1SignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the VAR1 model; deferred until port integration.
      enum class Var1SignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the VAR1 model; deferred until implementation.
      enum class Var1MonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for VAR1.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Var1Data : public ComponentData<real_type,
                                             index_type,
                                             Var1Parameters,
                                             Var1Buses,
                                             Var1SignalInputs,
                                             Var1SignalOutputs,
                                             Var1MonitorableVariables>
      {
        Var1Data() = default;

        using Parameters           = Var1Parameters;
        using Buses                = Var1Buses;
        using SignalInputs         = Var1SignalInputs;
        using SignalOutputs        = Var1SignalOutputs;
        using MonitorableVariables = Var1MonitorableVariables;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
