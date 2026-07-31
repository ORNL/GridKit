/**
 * @file Ac9cData.hpp
 * @brief Modeling data for the AC9C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC9C model.
      enum class Ac9cParameters
      {
        RC,      ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,      ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,      ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPR,     ///< \f$K_{\mathrm{PR}}\f$ Voltage regulator proportional gain
        KIR,     ///< \f$K_{\mathrm{IR}}\f$ Voltage regulator integral gain
        KDR,     ///< \f$K_{\mathrm{DR}}\f$ Voltage regulator derivative gain
        TDR,     ///< \f$T_{\mathrm{DR}}\f$ Lag time constant for derivative channel of PID controller
        VPIDmax, ///< \f$V_{\mathrm{PIDmax}}\f$ Maximum voltage regulator output
        VPIDmin, ///< \f$V_{\mathrm{PIDmin}}\f$ Minimum voltage regulator output
        KPA,     ///< \f$K_{\mathrm{PA}}\f$ Field current regulator proportional gain
        KIA,     ///< \f$K_{\mathrm{IA}}\f$ Field current regulator integral gain
        VAmax,   ///< \f$V_{\mathrm{Amax}}\f$ Maximum current regulator output
        VAmin,   ///< \f$V_{\mathrm{Amin}}\f$ Minimum current regulator output
        KA,      ///< \f$K_{\mathrm{A}}\f$ Controlled rectifier bridge equivalent gain
        TA,      ///< \f$T_{\mathrm{A}}\f$ Controlled rectifier bridge equivalent time constant
        VRmax,   ///< \f$V_{\mathrm{Rmax}}\f$ Maximum rectifier bridge output
        VRmin,   ///< \f$V_{\mathrm{Rmin}}\f$ Minimum rectifier bridge output
        KF,      ///< \f$K_{\mathrm{F}}\f$ Exciter field current feedback gain
        TF,      ///< \f$T_{\mathrm{F}}\f$ Field current feedback time constant
        KFW,     ///< \f$K_{\mathrm{FW}}\f$ Free wheel equivalent feedback gain
        VFWmax,  ///< \f$V_{\mathrm{FWmax}}\f$ Maximum free wheel feedback
        VFWmin,  ///< \f$V_{\mathrm{FWmin}}\f$ Minimum free wheel feedback
        SCT,     ///< \f$S_{\mathrm{CT}}\f$ Power stage type selector
        KC,      ///< \f$K_{\mathrm{C}}\f$ Diode bridge loading factor proportional to commutating reactance
        KD,      ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,      ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,      ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        VFEmax,  ///< \f$V_{\mathrm{FEmax}}\f$ Exciter field current limit
        VEmin,   ///< \f$V_{\mathrm{Emin}}\f$ Minimum exciter output limit
        E1,      ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,     ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,      ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2,     ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
        SW1,     ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KC1,     ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        KP,      ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI1,     ///< \f$K_{\mathrm{I1}}\f$ Compound circuit (current) gain coefficient
        XL,      ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP,  ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VB1max,  ///< \f$V_{\mathrm{B1max}}\f$ Maximum available exciter voltage
        KC2,     ///< \f$K_{\mathrm{C2}}\f$ Rectifier loading factor proportional to commutating reactance
        KI2,     ///< \f$K_{\mathrm{I2}}\f$ Compound circuit (current) gain coefficient
        VB2max   ///< \f$V_{\mathrm{B2max}}\f$ Maximum available compound exciter voltage
      };

      /// Bus keys for the AC9C model; deferred until port integration.
      enum class Ac9cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC9C model; deferred until port integration.
      enum class Ac9cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC9C model; deferred until port integration.
      enum class Ac9cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC9C model; deferred until implementation.
      enum class Ac9cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC9C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac9cData : public ComponentData<real_type,
                                             index_type,
                                             Ac9cParameters,
                                             Ac9cBuses,
                                             Ac9cSignalInputs,
                                             Ac9cSignalOutputs,
                                             Ac9cMonitorableVariables>
      {
        Ac9cData() = default;

        using Parameters           = Ac9cParameters;
        using Buses                = Ac9cBuses;
        using SignalInputs         = Ac9cSignalInputs;
        using SignalOutputs        = Ac9cSignalOutputs;
        using MonitorableVariables = Ac9cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
