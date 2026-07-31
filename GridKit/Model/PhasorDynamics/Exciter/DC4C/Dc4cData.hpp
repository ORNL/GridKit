/**
 * @file Dc4cData.hpp
 * @brief Modeling data for the DC4C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the DC4C model.
      enum class Dc4cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPR,    ///< \f$K_{\mathrm{PR}}\f$ Regulator proportional gain
        KIR,    ///< \f$K_{\mathrm{IR}}\f$ Regulator integral gain
        KDR,    ///< \f$K_{\mathrm{DR}}\f$ Regulator derivative gain
        TDR,    ///< \f$T_{\mathrm{DR}}\f$ Regulator derivative filter time constant
        KA,     ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum controller output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum controller output
        TE,     ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        KE,     ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        VEmin,  ///< \f$V_{\mathrm{Emin}}\f$ Exciter minimum output voltage
        E1,     ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        SE1,    ///< \f$S_{\mathrm{E1}}\f$ Exciter saturation factor at exciter output voltage E 1
        E2,     ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        SE2,    ///< \f$S_{\mathrm{E2}}\f$ Exciter saturation factor at exciter output voltage E 2
        KF,     ///< \f$K_{\mathrm{F}}\f$ Rate feedback gain
        TF,     ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit gain coefficient
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        KI,     ///< \f$K_{\mathrm{I}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        KC1,    ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        VBmax,  ///< \f$V_{\mathrm{Bmax}}\f$ Maximum available exciter field voltage
        SW1     ///< \f$\mathrm{SW}_{1}\f$ Logical switch 1
      };

      /// Bus keys for the DC4C model; deferred until port integration.
      enum class Dc4cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the DC4C model; deferred until port integration.
      enum class Dc4cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the DC4C model; deferred until port integration.
      enum class Dc4cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the DC4C model; deferred until implementation.
      enum class Dc4cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for DC4C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Dc4cData : public ComponentData<real_type,
                                             index_type,
                                             Dc4cParameters,
                                             Dc4cBuses,
                                             Dc4cSignalInputs,
                                             Dc4cSignalOutputs,
                                             Dc4cMonitorableVariables>
      {
        Dc4cData() = default;

        using Parameters           = Dc4cParameters;
        using Buses                = Dc4cBuses;
        using SignalInputs         = Dc4cSignalInputs;
        using SignalOutputs        = Dc4cSignalOutputs;
        using MonitorableVariables = Dc4cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
