/**
 * @file St7cData.hpp
 * @brief Modeling data for the ST7C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the ST7C model.
      enum class St7cParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        TG,    ///< \f$T_{\mathrm{G}}\f$ Regulator input filter time constant
        TF,    ///< \f$T_{\mathrm{F}}\f$ Voltage regulator denominator (lag) time constant
        Vmax,  ///< \f$V_{\mathrm{max}}\f$ Maximum voltage reference
        Vmin,  ///< \f$V_{\mathrm{min}}\f$ Minimum voltage reference
        KPA,   ///< \f$K_{\mathrm{PA}}\f$ Voltage regulator gain
        TB,    ///< \f$T_{\mathrm{B}}\f$ Voltage regulator denominator (lag) time constant
        TC,    ///< \f$T_{\mathrm{C}}\f$ Voltage regulator numerator (lead) time constant
        TA,    ///< \f$T_{\mathrm{A}}\f$ Thyristor bridge firing control equivalent time constant
        VRmax, ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin, ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KL,    ///< \f$K_{\mathrm{L}}\f$ Minimum excitation limit gain
        KH,    ///< \f$K_{\mathrm{H}}\f$ Maximum excitation limit gain
        KIA,   ///< \f$K_{\mathrm{IA}}\f$ PI regulator feedback gain
        TIA    ///< \f$T_{\mathrm{IA}}\f$ PI regulator feedback time constant
      };

      /// Bus keys for the ST7C model; deferred until port integration.
      enum class St7cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the ST7C model; deferred until port integration.
      enum class St7cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the ST7C model; deferred until port integration.
      enum class St7cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the ST7C model; deferred until implementation.
      enum class St7cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for ST7C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct St7cData : public ComponentData<real_type,
                                             index_type,
                                             St7cParameters,
                                             St7cBuses,
                                             St7cSignalInputs,
                                             St7cSignalOutputs,
                                             St7cMonitorableVariables>
      {
        St7cData() = default;

        using Parameters           = St7cParameters;
        using Buses                = St7cBuses;
        using SignalInputs         = St7cSignalInputs;
        using SignalOutputs        = St7cSignalOutputs;
        using MonitorableVariables = St7cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
