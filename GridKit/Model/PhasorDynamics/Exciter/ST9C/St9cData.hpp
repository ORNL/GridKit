/**
 * @file St9cData.hpp
 * @brief Modeling data for the ST9C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST9C model.
      enum class St9cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        TCD,    ///< \f$T_{\mathrm{CD}}\f$ Time constant of differential part of AVR
        TBD,    ///< \f$T_{\mathrm{BD}}\f$ Filter time constant of differential part of AVR
        ZA,     ///< \f$Z_{\mathrm{A}}\f$ Dead-band for differential part influence on AVR
        KA,     ///< \f$K_{\mathrm{A}}\f$ AVR gain
        KU,     ///< \f$K_{\mathrm{U}}\f$ Gain associated with activation of takeover UEL
        TA,     ///< \f$T_{\mathrm{A}}\f$ Time constant of AVR
        TAUEL,  ///< \f$T_{\mathrm{AUEL}}\f$ Time constant of underexcitation limiter
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        KAS,    ///< \f$K_{\mathrm{AS}}\f$ Power converter gain (proportional to supply voltage)
        TAS,    ///< \f$T_{\mathrm{AS}}\f$ Equivalent time constant of power converter firing control
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with the compound source
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Switch to select exciter supply configuration
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        VBmax   ///< \f$V_{\mathrm{Bmax}}\f$ Maximum limit on exciter voltage based on supply condition
      };

      /// Bus keys for the ST9C model; deferred until port integration.
      enum class St9cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST9C model; deferred until port integration.
      enum class St9cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST9C model; deferred until port integration.
      enum class St9cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST9C model; deferred until implementation.
      enum class St9cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST9C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St9cData : public ComponentData<real_type,
                                             index_type,
                                             St9cParameters,
                                             St9cBuses,
                                             St9cSignalInputs,
                                             St9cSignalOutputs,
                                             St9cMonitorableVariables>
      {
        St9cData() = default;

        using Parameters           = St9cParameters;
        using Buses                = St9cBuses;
        using SignalInputs         = St9cSignalInputs;
        using SignalOutputs        = St9cSignalOutputs;
        using MonitorableVariables = St9cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
