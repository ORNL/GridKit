/**
 * @file St1cData.hpp
 * @brief Modeling data for the ST1C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST1C model.
      enum class St1cParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,    ///< \f$K_{\mathrm{A}}\f$ Voltage regulator gain
        TA,    ///< \f$T_{\mathrm{A}}\f$ Voltage regulator time constant
        TB,    ///< \f$T_{\mathrm{B}}\f$ Regulator denominator (lag) time constant
        TC,    ///< \f$T_{\mathrm{C}}\f$ Regulator numerator (lead) time constant
        TB1,   ///< \f$T_{\mathrm{B1}}\f$ Regulator denominator (lag) time constant
        TC1,   ///< \f$T_{\mathrm{C1}}\f$ Regulator numerator (lead) time constant
        VRmax, ///< \f$V_{\mathrm{Rmax}}\f$ Maximum exciter output
        VRmin, ///< \f$V_{\mathrm{Rmin}}\f$ Minimum exciter output
        VAmax, ///< \f$V_{\mathrm{Amax}}\f$ Maximum regulator output
        VAmin, ///< \f$V_{\mathrm{Amin}}\f$ Minimum regulator output
        VImax, ///< \f$V_{\mathrm{Imax}}\f$ Maximum voltage error (regulator input)
        VImin, ///< \f$V_{\mathrm{Imin}}\f$ Minimum voltage error (regulator input)
        KF,    ///< \f$K_{\mathrm{F}}\f$ Rate feedback gain
        TF,    ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        KC,    ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KLR,   ///< \f$K_{\mathrm{LR}}\f$ Exciter output current limiter gain
        ILR    ///< \f$I_{\mathrm{LR}}\f$ Exciter output current limit reference
      };

      /// Bus keys for the ST1C model; deferred until port integration.
      enum class St1cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST1C model; deferred until port integration.
      enum class St1cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST1C model; deferred until port integration.
      enum class St1cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST1C model; deferred until implementation.
      enum class St1cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST1C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St1cData : public ComponentData<real_type,
                                             index_type,
                                             St1cParameters,
                                             St1cBuses,
                                             St1cSignalInputs,
                                             St1cSignalOutputs,
                                             St1cMonitorableVariables>
      {
        St1cData() = default;

        using Parameters           = St1cParameters;
        using Buses                = St1cBuses;
        using SignalInputs         = St1cSignalInputs;
        using SignalOutputs        = St1cSignalOutputs;
        using MonitorableVariables = St1cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
