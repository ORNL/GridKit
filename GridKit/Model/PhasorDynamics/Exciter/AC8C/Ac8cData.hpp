/**
 * @file Ac8cData.hpp
 * @brief Modeling data for the AC8C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC8C model.
      enum class Ac8cParameters
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
        KA,     ///< \f$K_{\mathrm{A}}\f$ Rectifier bridge gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Rectifier bridge time constant
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit gain coefficient
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        KC1,    ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        VBmax,  ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter field voltage
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Logical switch 1
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KD,     ///< \f$K_{\mathrm{D}}\f$ Demagnetizing factor, function of exciter alternator reactances
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        VFEMAX, ///< \f$V_{\mathrm{FEMAX}}\f$ Maximum exciter field current
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,    ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2     ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the AC8C model; deferred until port integration.
      enum class Ac8cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC8C model; deferred until port integration.
      enum class Ac8cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC8C model; deferred until port integration.
      enum class Ac8cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC8C model; deferred until implementation.
      enum class Ac8cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC8C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac8cData : public ComponentData<real_type,
                                             index_type,
                                             Ac8cParameters,
                                             Ac8cBuses,
                                             Ac8cSignalInputs,
                                             Ac8cSignalOutputs,
                                             Ac8cMonitorableVariables>
      {
        Ac8cData() = default;

        using Parameters           = Ac8cParameters;
        using Buses                = Ac8cBuses;
        using SignalInputs         = Ac8cSignalInputs;
        using SignalOutputs        = Ac8cSignalOutputs;
        using MonitorableVariables = Ac8cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
