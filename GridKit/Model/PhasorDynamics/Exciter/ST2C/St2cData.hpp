/**
 * @file St2cData.hpp
 * @brief Modeling data for the ST2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST2C model.
      enum class St2cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,     ///< \f$K_{\mathrm{A}}\f$ Voltage regulator gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Voltage regulator time constant
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KF,     ///< \f$K_{\mathrm{F}}\f$ Rate feedback gain
        TF,     ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        KC,     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        EFDmax, ///< \f$E_{\mathrm{FDmax}}\f$ Maximum generator field voltage
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI,     ///< \f$K_{\mathrm{I}}\f$ Compound circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VBmax   ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the ST2C model; deferred until port integration.
      enum class St2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST2C model; deferred until port integration.
      enum class St2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST2C model; deferred until port integration.
      enum class St2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST2C model; deferred until implementation.
      enum class St2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St2cData : public ComponentData<real_type,
                                             index_type,
                                             St2cParameters,
                                             St2cBuses,
                                             St2cSignalInputs,
                                             St2cSignalOutputs,
                                             St2cMonitorableVariables>
      {
        St2cData() = default;

        using Parameters           = St2cParameters;
        using Buses                = St2cBuses;
        using SignalInputs         = St2cSignalInputs;
        using SignalOutputs        = St2cSignalOutputs;
        using MonitorableVariables = St2cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
