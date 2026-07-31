/**
 * @file St5cData.hpp
 * @brief Modeling data for the ST5C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST5C model.
      enum class St5cParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KR,    ///< \f$K_{\mathrm{R}}\f$ Regulator gain
        TB1,   ///< \f$T_{\mathrm{B1}}\f$ Voltage regulator denominator (lag) time constant (first block)
        TC1,   ///< \f$T_{\mathrm{C1}}\f$ Voltage regulator numerator (lead) time constant (first block)
        TB2,   ///< \f$T_{\mathrm{B2}}\f$ Voltage regulator denominator (lag) time constant (second block)
        TC2,   ///< \f$T_{\mathrm{C2}}\f$ Voltage regulator numerator (lead) time constant (second block)
        TUB1,  ///< \f$T_{\mathrm{UB1}}\f$ UEL regulator denominator (lag) time constant (first block)
        TUC1,  ///< \f$T_{\mathrm{UC1}}\f$ UEL regulator numerator (lead) time constant (first block)
        TUB2,  ///< \f$T_{\mathrm{UB2}}\f$ UEL regulator denominator (lag) time constant (second block)
        TUC2,  ///< \f$T_{\mathrm{UC2}}\f$ UEL regulator numerator (lead) time constant (second block)
        TOB1,  ///< \f$T_{\mathrm{OB1}}\f$ OEL regulator denominator (lag) time constant (first block)
        TOC1,  ///< \f$T_{\mathrm{OC1}}\f$ OEL regulator numerator (lead) time constant (first block)
        TOB2,  ///< \f$T_{\mathrm{OB2}}\f$ OEL regulator denominator (lag) time constant (second block)
        TOC2,  ///< \f$T_{\mathrm{OC2}}\f$ OEL regulator numerator (lead) time constant (second block)
        VRmax, ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin, ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KC,    ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        T1     ///< \f$T_{1}\f$ Thyristor bridge firing control equivalent time constant
      };

      /// Bus keys for the ST5C model; deferred until port integration.
      enum class St5cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST5C model; deferred until port integration.
      enum class St5cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST5C model; deferred until port integration.
      enum class St5cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST5C model; deferred until implementation.
      enum class St5cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST5C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St5cData : public ComponentData<real_type,
                                             index_type,
                                             St5cParameters,
                                             St5cBuses,
                                             St5cSignalInputs,
                                             St5cSignalOutputs,
                                             St5cMonitorableVariables>
      {
        St5cData() = default;

        using Parameters           = St5cParameters;
        using Buses                = St5cBuses;
        using SignalInputs         = St5cSignalInputs;
        using SignalOutputs        = St5cSignalOutputs;
        using MonitorableVariables = St5cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
