/**
 * @file Ac10cData.hpp
 * @brief Modeling data for the AC10C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC10C model.
      enum class Ac10cParameters
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
        SWEXC,  ///< \f$\mathrm{SW}_{\mathrm{EXC}}\f$ Exciter field current regulator feedback selector
        TEXC,   ///< \f$T_{\mathrm{EXC}}\f$ Exciter field current regulator measurement time constant
        KEXC,   ///< \f$K_{\mathrm{EXC}}\f$ Exciter field current regulator feedback gain
        KCR,    ///< \f$K_{\mathrm{CR}}\f$ Exciter field current regulator proportional gain
        TF1,    ///< \f$T_{\mathrm{F1}}\f$ Exciter field current regulator numerator (lead) time constant
        TF2,    ///< \f$T_{\mathrm{F2}}\f$ Exciter field current regulator denominator (lag) time constant
        KVFE,   ///< \f$K_{\mathrm{VFE}}\f$ Exciter field current limiter feedback gain
        KLIM,   ///< \f$K_{\mathrm{LIM}}\f$ Exciter field current limiter proportional gain
        VFELIM, ///< \f$V_{\mathrm{FELIM}}\f$ Exciter field current limiter reference
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KD,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        VFEmax, ///< \f$V_{\mathrm{FEmax}}\f$ Maximum exciter field current
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E) E 1
        Se1,    ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E) E 2
        Se2,    ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        KC1,    ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        VB1max, ///< \f$V_{\mathrm{B1max}}\f$ Maximum available exciter voltage
        KI2,    ///< \f$K_{\mathrm{I2}}\f$ Additive potential circuit (current) gain coefficient
        KC2,    ///< \f$K_{\mathrm{C2}}\f$ Rectifier loading factor proportional to commutating reactance
        VB2max  ///< \f$V_{\mathrm{B2max}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the AC10C model; deferred until port integration.
      enum class Ac10cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC10C model; deferred until port integration.
      enum class Ac10cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC10C model; deferred until port integration.
      enum class Ac10cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC10C model; deferred until implementation.
      enum class Ac10cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC10C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac10cData : public ComponentData<real_type,
                                              index_type,
                                              Ac10cParameters,
                                              Ac10cBuses,
                                              Ac10cSignalInputs,
                                              Ac10cSignalOutputs,
                                              Ac10cMonitorableVariables>
      {
        Ac10cData() = default;

        using Parameters           = Ac10cParameters;
        using Buses                = Ac10cBuses;
        using SignalInputs         = Ac10cSignalInputs;
        using SignalOutputs        = Ac10cSignalOutputs;
        using MonitorableVariables = Ac10cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
