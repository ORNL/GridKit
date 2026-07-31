/**
 * @file St4cData.hpp
 * @brief Modeling data for the ST4C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST4C model.
      enum class St4cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPR,    ///< \f$K_{\mathrm{PR}}\f$ Voltage regulator proportional gain
        KIR,    ///< \f$K_{\mathrm{IR}}\f$ Voltage regulator integral gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Thyristor bridge firing control equivalent time constant
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KPM,    ///< \f$K_{\mathrm{PM}}\f$ Forward proportional gain of inner loop field regulator
        KIM,    ///< \f$K_{\mathrm{IM}}\f$ Forward integral gain of inner loop field regulator
        VMmax,  ///< \f$V_{\mathrm{Mmax}}\f$ Maximum output of inner loop field regulator
        VMmin,  ///< \f$V_{\mathrm{Mmin}}\f$ Minimum output of inner loop field regulator
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum exciter output
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum exciter output
        KG,     ///< \f$K_{\mathrm{G}}\f$ Feedback gain of inner loop field regulator
        TG,     ///< \f$T_{\mathrm{G}}\f$ Feedback time constant of field current regulator
        VGmax,  ///< \f$V_{\mathrm{Gmax}}\f$ Maximum feedback voltage for field current regulator
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential or compound circuit current gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VBmax   ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the ST4C model; deferred until port integration.
      enum class St4cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST4C model; deferred until port integration.
      enum class St4cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST4C model; deferred until port integration.
      enum class St4cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST4C model; deferred until implementation.
      enum class St4cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST4C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St4cData : public ComponentData<real_type,
                                             index_type,
                                             St4cParameters,
                                             St4cBuses,
                                             St4cSignalInputs,
                                             St4cSignalOutputs,
                                             St4cMonitorableVariables>
      {
        St4cData() = default;

        using Parameters           = St4cParameters;
        using Buses                = St4cBuses;
        using SignalInputs         = St4cSignalInputs;
        using SignalOutputs        = St4cSignalOutputs;
        using MonitorableVariables = St4cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
