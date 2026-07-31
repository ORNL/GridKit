/**
 * @file St10cData.hpp
 * @brief Modeling data for the ST10C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST10C model.
      enum class St10cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KR,     ///< \f$K_{\mathrm{R}}\f$ Regulator gain
        TB1,    ///< \f$T_{\mathrm{B1}}\f$ Voltage regulator denominator (lag) time constant 1
        TC1,    ///< \f$T_{\mathrm{C1}}\f$ Voltage regulator numerator (lead) time constant 1
        TB2,    ///< \f$T_{\mathrm{B2}}\f$ Voltage regulator denominator (lag) time constant 2
        TC2,    ///< \f$T_{\mathrm{C2}}\f$ Voltage regulator numerator (lead) time constant 2
        TUB1,   ///< \f$T_{\mathrm{UB1}}\f$ UEL regulator denominator (lag) time constant 1
        TUC1,   ///< \f$T_{\mathrm{UC1}}\f$ UEL regulator numerator (lead) time constant 1
        TUB2,   ///< \f$T_{\mathrm{UB2}}\f$ UEL regulator denominator (lag) time constant 2
        TUC2,   ///< \f$T_{\mathrm{UC2}}\f$ UEL regulator numerator (lead) time constant 2
        TOB1,   ///< \f$T_{\mathrm{OB1}}\f$ OEL regulator denominator (lag) time constant 1
        TOC1,   ///< \f$T_{\mathrm{OC1}}\f$ OEL regulator numerator (lead) time constant 1
        TOB2,   ///< \f$T_{\mathrm{OB2}}\f$ OEL regulator denominator (lag) time constant 2
        TOC2,   ///< \f$T_{\mathrm{OC2}}\f$ OEL regulator numerator (lead) time constant 2
        VRSmax, ///< \f$V_{\mathrm{RSmax}}\f$ Maximum PSS regulator output
        VRSmin, ///< \f$V_{\mathrm{RSmin}}\f$ Minimum PSS regulator output
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        T1,     ///< \f$T_{1}\f$ Equivalent time constant for rectifier bridge
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VBmax   ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the ST10C model; deferred until port integration.
      enum class St10cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST10C model; deferred until port integration.
      enum class St10cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST10C model; deferred until port integration.
      enum class St10cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST10C model; deferred until implementation.
      enum class St10cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST10C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St10cData : public ComponentData<real_type,
                                              index_type,
                                              St10cParameters,
                                              St10cBuses,
                                              St10cSignalInputs,
                                              St10cSignalOutputs,
                                              St10cMonitorableVariables>
      {
        St10cData() = default;

        using Parameters           = St10cParameters;
        using Buses                = St10cBuses;
        using SignalInputs         = St10cSignalInputs;
        using SignalOutputs        = St10cSignalOutputs;
        using MonitorableVariables = St10cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
