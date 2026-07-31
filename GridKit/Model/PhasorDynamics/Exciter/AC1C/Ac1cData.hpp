/**
 * @file Ac1cData.hpp
 * @brief Modeling data for the AC1C exciter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC1C exciter model.
      enum class Ac1cParameters
      {
        Rc,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        Xc,     ///< \f$X_{\mathrm{C}}\f$ Reactive component of load compensation
        Tr,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        Ka,     ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        Ta,     ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        Tb,     ///< \f$T_{\mathrm{B}}\f$ Regulator denominator time constant
        Tc,     ///< \f$T_{\mathrm{C}}\f$ Regulator numerator time constant
        Kf,     ///< \f$K_{\mathrm{F}}\f$ Rate-feedback stabilizer gain
        Tf,     ///< \f$T_{\mathrm{F}}\f$ Rate-feedback time constant
        Vamax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum regulator output
        Vamin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum regulator output
        Efemax, ///< \f$E_{\mathrm{FEmax}}\f$ Maximum exciter field voltage
        Efemin, ///< \f$E_{\mathrm{FEmin}}\f$ Minimum exciter field voltage
        Ke,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        Te,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        Kc,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor
        Kd,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor
        E1,     ///< \f$E_1\f$ First exciter saturation reference voltage
        Se1,    ///< \f$S_{\mathrm{E}}(E_1)\f$ Saturation at \f$E_1\f$
        E2,     ///< \f$E_2\f$ Second exciter saturation reference voltage
        Se2     ///< \f$S_{\mathrm{E}}(E_2)\f$ Saturation at \f$E_2\f$
      };

      /// Buses for the AC1C exciter model.
      enum class Ac1cBuses : size_t
      {
        bus, ///< Required terminal-bus ID
        SIZE
      };

      /// Signal inputs for the AC1C exciter model.
      enum class Ac1cSignalInputs : size_t
      {
        ifd,  ///< Required generator field-current signal ID
        ir,   ///< Optional terminal-current real-part signal ID for load compensation
        ii,   ///< Optional terminal-current imaginary-part signal ID for load compensation
        vref, ///< Optional voltage-reference signal ID
        vs,   ///< Optional stabilizer output signal ID
        vuel, ///< Optional underexcitation-limiter output signal ID
        voel, ///< Optional overexcitation-limiter output signal ID
        vscl, ///< Optional stator-current-limiter output signal ID
        SIZE
      };

      /// Signal outputs for the AC1C exciter model.
      enum class Ac1cSignalOutputs : size_t
      {
        efd, ///< Required generator field-voltage signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class Ac1cMonitorableVariables
      {
        efd, ///< Generator field-voltage output
        ve,  ///< Exciter voltage
        vts, ///< Sensed terminal voltage
        va,  ///< Voltage-regulator output
        vll, ///< Lead-lag block output
        vf   ///< Rate-feedback signal
      };

      /**
       * @brief Model data for AC1C: parameters, the terminal bus, signal
       *        connections, and monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac1cData : public ComponentData<real_type,
                                             index_type,
                                             Ac1cParameters,
                                             Ac1cBuses,
                                             Ac1cSignalInputs,
                                             Ac1cSignalOutputs,
                                             Ac1cMonitorableVariables>
      {
        Ac1cData() = default;

        using Parameters           = Ac1cParameters;
        using Buses                = Ac1cBuses;
        using SignalInputs         = Ac1cSignalInputs;
        using SignalOutputs        = Ac1cSignalOutputs;
        using MonitorableVariables = Ac1cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
