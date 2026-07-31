/**
 * @file St8cData.hpp
 * @brief Modeling data for the ST8C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST8C model.
      enum class St8cParameters
      {
        RC,     ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,     ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,     ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KPR,    ///< \f$K_{\mathrm{PR}}\f$ Voltage regulator proportional gain
        KIR,    ///< \f$K_{\mathrm{IR}}\f$ Voltage regulator integral gain
        VPImax, ///< \f$V_{\mathrm{PImax}}\f$ Maximum voltage regulator output
        VPImin, ///< \f$V_{\mathrm{PImin}}\f$ Minimum voltage regulator output
        KPA,    ///< \f$K_{\mathrm{PA}}\f$ Field current regulator proportional gain
        KIA,    ///< \f$K_{\mathrm{IA}}\f$ Field current regulator integral gain
        VAmax,  ///< \f$V_{\mathrm{Amax}}\f$ Maximum field current regulator output
        VAmin,  ///< \f$V_{\mathrm{Amin}}\f$ Minimum field current regulator output
        KA,     ///< \f$K_{\mathrm{A}}\f$ Field current regulator proportional gain
        TA,     ///< \f$T_{\mathrm{A}}\f$ Controlled rectifier bridge equivalent time constant
        VRmax,  ///< \f$V_{\mathrm{Rmax}}\f$ Maximum field current regulator output
        VRmin,  ///< \f$V_{\mathrm{Rmin}}\f$ Minimum field current regulator output
        KF,     ///< \f$K_{\mathrm{F}}\f$ Exciter field current feedback gain
        TF,     ///< \f$T_{\mathrm{F}}\f$ Field current feedback time constant
        SW1,    ///< \f$\mathrm{SW}_{1}\f$ Power source selector
        KC1,    ///< \f$K_{\mathrm{C1}}\f$ Rectifier loading factor proportional to commutating reactance
        KP,     ///< \f$K_{\mathrm{P}}\f$ Potential circuit (voltage) gain coefficient
        KI1,    ///< \f$K_{\mathrm{I1}}\f$ Potential circuit (current) gain coefficient
        XL,     ///< \f$X_{\mathrm{L}}\f$ Reactance associated with potential source
        ThetaP, ///< \f$\theta_{\mathrm{P}}\f$ Potential circuit phase angle (degrees)
        VBmax1, ///< \f$V_{\mathrm{Bmax1}}\f$ Maximum available exciter voltage
        KC2,    ///< \f$K_{\mathrm{C2}}\f$ Rectifier loading factor proportional to commutating reactance
        KI2,    ///< \f$K_{\mathrm{I2}}\f$ Potential circuit (current) gain coefficient
        VBmax2  ///< \f$V_{\mathrm{Bmax2}}\f$ Maximum available exciter voltage
      };

      /// Bus keys for the ST8C model; deferred until port integration.
      enum class St8cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST8C model; deferred until port integration.
      enum class St8cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST8C model; deferred until port integration.
      enum class St8cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST8C model; deferred until implementation.
      enum class St8cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST8C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St8cData : public ComponentData<real_type,
                                             index_type,
                                             St8cParameters,
                                             St8cBuses,
                                             St8cSignalInputs,
                                             St8cSignalOutputs,
                                             St8cMonitorableVariables>
      {
        St8cData() = default;

        using Parameters           = St8cParameters;
        using Buses                = St8cBuses;
        using SignalInputs         = St8cSignalInputs;
        using SignalOutputs        = St8cSignalOutputs;
        using MonitorableVariables = St8cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
