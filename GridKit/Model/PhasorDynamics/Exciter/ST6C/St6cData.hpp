/**
 * @file St6cData.hpp
 * @brief Modeling data for the ST6C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST6C model.
      enum class St6cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPA,    ///< \f$K_{\mathrm{PA}}\f$ Voltage regulator proportional gain
        KIA,    ///< \f$K_{\mathrm{IA}}\f$ Voltage regulator integral gain
        KFF,    ///< \f$K_{\mathrm{FF}}\f$ Pre-control gain constant of the inner loop field voltage regulator
        KM,     ///< \f$K_{\mathrm{M}}\f$ Forward gain constant of the inner loop field voltage regulator
        KG,     ///< \f$K_{\mathrm{G}}\f$ Feedback gain constant of the inner loop field voltage regulator
        TG,     ///< \f$T_{\mathrm{G}}\f$ Feedback time constant of inner loop field voltage regulator
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum voltage regulator output
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum voltage regulator output
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output limit
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output limit
        VMmax,  ///< \f$V_{\mathrm{Mmax}}\f$ Maximum rectifier output limit
        VMmin,  ///< \f$V_{\mathrm{Mmin}}\f$ Minimum rectifier output limit
        TA,     ///< \f$T_{\mathrm{A}}\f$ Thyristor bridge firing control equivalent time constant
        KLR,    ///< \f$K_{\mathrm{LR}}\f$ Exciter output current limiter gain
        KCI,    ///< \f$K_{\mathrm{CI}}\f$ Exciter output current limit adjustment
        ILR,    ///< \f$I_{\mathrm{LR}}\f$ Exciter output current limit reference
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VBmax   ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the ST6C model; deferred until port integration.
      enum class St6cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST6C model; deferred until port integration.
      enum class St6cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST6C model; deferred until port integration.
      enum class St6cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST6C model; deferred until implementation.
      enum class St6cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST6C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St6cData : public ComponentData<real_type,
                                             index_type,
                                             St6cParameters,
                                             St6cBuses,
                                             St6cSignalInputs,
                                             St6cSignalOutputs,
                                             St6cMonitorableVariables>
      {
        St6cData() = default;

        using Parameters           = St6cParameters;
        using Buses                = St6cBuses;
        using SignalInputs         = St6cSignalInputs;
        using SignalOutputs        = St6cSignalOutputs;
        using MonitorableVariables = St6cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
