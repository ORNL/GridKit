/**
 * @file Ac7cData.hpp
 * @brief Modeling data for the AC7C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC7C model.
      enum class Ac7cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPR,    ///< \f$K_{\mathrm{PR}}\f$ Voltage regulator proportional gain
        KIR,    ///< \f$K_{\mathrm{IR}}\f$ Voltage regulator integral gain
        KDR,    ///< \f$K_{\mathrm{DR}}\f$ Voltage regulator derivative gain
        TDR,    ///< \f$T_{\mathrm{DR}}\f$ Lag time constant for derivative channel of PID controller
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KPA,    ///< \f$K_{\mathrm{PA}}\f$ Field current regulator proportional gain
        KIA,    ///< \f$K_{\mathrm{IA}}\f$ Field current regulator integral gain
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum field current regulator output
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum field current regulator output
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit gain coefficient
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        KC1,    ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        VBmax,  ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter field voltage
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        SW2,    ///< \f$\mathrm{SW}_{2}\f$ Power source selector
        KR,     ///< \f$K_{\mathrm{R}}\f$ Gain related to regulator and alternator field power supply
        KL,     ///< \f$K_{\mathrm{L}}\f$ Gain related to negative exciter field current capability
        KF1,    ///< \f$K_{\mathrm{F1}}\f$ Generator field voltage feedback gain
        KF2,    ///< \f$K_{\mathrm{F2}}\f$ Exciter field current feedback gain
        KF3,    ///< \f$K_{\mathrm{F3}}\f$ Rate feedback gain
        TF,     ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KD,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        VFEmax, ///< \f$V_{\mathrm{FEmax}}\f$ Maximum exciter field current
        VEmin,  ///< \f$V_{\mathrm{Emin}}\f$ Minimum exciter voltage output
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,    ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2     ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the AC7C model; deferred until port integration.
      enum class Ac7cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC7C model; deferred until port integration.
      enum class Ac7cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC7C model; deferred until port integration.
      enum class Ac7cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC7C model; deferred until implementation.
      enum class Ac7cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC7C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac7cData : public ComponentData<real_type,
                                             index_type,
                                             Ac7cParameters,
                                             Ac7cBuses,
                                             Ac7cSignalInputs,
                                             Ac7cSignalOutputs,
                                             Ac7cMonitorableVariables>
      {
        Ac7cData() = default;

        using Parameters           = Ac7cParameters;
        using Buses                = Ac7cBuses;
        using SignalInputs         = Ac7cSignalInputs;
        using SignalOutputs        = Ac7cSignalOutputs;
        using MonitorableVariables = Ac7cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
