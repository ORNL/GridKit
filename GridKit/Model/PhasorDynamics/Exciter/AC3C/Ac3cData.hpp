/**
 * @file Ac3cData.hpp
 * @brief Modeling data for the AC3C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC3C model.
      enum class Ac3cParameters
      {
        RC,      ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,      ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,      ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPR,     ///< \f$K_{\mathrm{PR}}\f$ Voltage regulator proportional gain
        KIR,     ///< \f$K_{\mathrm{IR}}\f$ Voltage regulator integral gain
        KDR,     ///< \f$K_{\mathrm{DR}}\f$ Voltage regulator derivative gain
        TDR,     ///< \f$T_{\mathrm{DR}}\f$ Lag time constant for derivative channel of PID controller
        VPIDmax, ///< \f$V_{\mathrm{PIDmax}}\f$ Maximum PID regulator output
        VPIDmin, ///< \f$V_{\mathrm{PIDmin}}\f$ Minimum PID regulator output
        TC,      ///< \f$T_{\mathrm{C}}\f$ Regulator numerator (lead) time constant
        TB,      ///< \f$T_{\mathrm{B}}\f$ Regulator denominator (lag) time constant
        KA,      ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        TA,      ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        VAmax,   ///< \f$V_{\mathrm{Amax}}\f$ Maximum regulator output
        VAmin,   ///< \f$V_{\mathrm{Amin}}\f$ Minimum regulator output
        TE,      ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        TF,      ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        EFDN,    ///< \f$E_{\mathrm{FDN}}\f$ Value of E at which feedback gain changes FD
        VEmin,   ///< \f$V_{\mathrm{Emin}}\f$ Minimum exciter voltage output
        E1,      ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,     ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,      ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2,     ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
        KR,      ///< \f$K_{\mathrm{R}}\f$ Gain associated with regulator and alternator field power supply
        KC,      ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KD,      ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,      ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        KF,      ///< \f$K_{\mathrm{F}}\f$ Rate feedback excitation system stabilizer gain
        KN,      ///< \f$K_{\mathrm{N}}\f$ Rate feedback excitation system stabilizer gain
        VFEmax   ///< \f$V_{\mathrm{FEmax}}\f$ Exciter field current limit
      };

      /// Bus keys for the AC3C model; deferred until port integration.
      enum class Ac3cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC3C model; deferred until port integration.
      enum class Ac3cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC3C model; deferred until port integration.
      enum class Ac3cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC3C model; deferred until implementation.
      enum class Ac3cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC3C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac3cData : public ComponentData<real_type,
                                             index_type,
                                             Ac3cParameters,
                                             Ac3cBuses,
                                             Ac3cSignalInputs,
                                             Ac3cSignalOutputs,
                                             Ac3cMonitorableVariables>
      {
        Ac3cData() = default;

        using Parameters           = Ac3cParameters;
        using Buses                = Ac3cBuses;
        using SignalInputs         = Ac3cSignalInputs;
        using SignalOutputs        = Ac3cSignalOutputs;
        using MonitorableVariables = Ac3cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
