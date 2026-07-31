/**
 * @file St3cData.hpp
 * @brief Modeling data for the ST3C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST3C model.
      enum class St3cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,     ///< \f$K_{\mathrm{A}}\f$ Voltage regulator gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Thyristor bridge firing control equivalent time constant
        TB,     ///< \f$T_{\mathrm{B}}\f$ Voltage regulator denominator (lag) time constant
        TC,     ///< \f$T_{\mathrm{C}}\f$ Voltage regulator numerator (lead) time constant
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum voltage regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum voltage regulator output
        VImax,  ///< \f$V_{\mathrm{Imax}}\f$ Maximum voltage error (regulator input)
        VImin,  ///< \f$V_{\mathrm{Imin}}\f$ Minimum voltage error (regulator input)
        KM,     ///< \f$K_{\mathrm{M}}\f$ Forward gain of inner loop field regulator
        TM,     ///< \f$T_{\mathrm{M}}\f$ Forward time constant of inner loop field regulator
        VMmax,  ///< \f$V_{\mathrm{Mmax}}\f$ Maximum output of field current regulator
        VMmin,  ///< \f$V_{\mathrm{Mmin}}\f$ Minimum output of field current regulator
        KG,     ///< \f$K_{\mathrm{G}}\f$ Feedback gain of field current regulator
        VGmax,  ///< \f$V_{\mathrm{Gmax}}\f$ Maximum feedback voltage for field current regulator
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Compound circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VBmax   ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the ST3C model; deferred until port integration.
      enum class St3cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST3C model; deferred until port integration.
      enum class St3cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST3C model; deferred until port integration.
      enum class St3cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST3C model; deferred until implementation.
      enum class St3cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST3C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St3cData : public ComponentData<real_type,
                                             index_type,
                                             St3cParameters,
                                             St3cBuses,
                                             St3cSignalInputs,
                                             St3cSignalOutputs,
                                             St3cMonitorableVariables>
      {
        St3cData() = default;

        using Parameters           = St3cParameters;
        using Buses                = St3cBuses;
        using SignalInputs         = St3cSignalInputs;
        using SignalOutputs        = St3cSignalOutputs;
        using MonitorableVariables = St3cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
