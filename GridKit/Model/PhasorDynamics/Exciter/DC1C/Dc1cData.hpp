/**
 * @file Dc1cData.hpp
 * @brief Modeling data for the DC1C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the DC1C model.
      enum class Dc1cParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,    ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        TA,    ///< \f$T_{\mathrm{A}}\f$ Regulator time constant
        TB,    ///< \f$T_{\mathrm{B}}\f$ Regulator denominator (lag) time constant
        TC,    ///< \f$T_{\mathrm{C}}\f$ Regulator numerator (lead) time constant
        KE,    ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        TE,    ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        VRMAX, ///< \f$V_{\mathrm{RMAX}}\f$ Maximum controller output
        VRMIN, ///< \f$V_{\mathrm{RMIN}}\f$ Minimum controller output
        KF,    ///< \f$K_{\mathrm{F}}\f$ Rate feedback gain
        TF,    ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        E1,    ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        SE1,   ///< \f$S_{\mathrm{E1}}\f$ Exciter saturation factor at exciter output voltage E 1
        E2,    ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        SE2    ///< \f$S_{\mathrm{E2}}\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the DC1C model; deferred until port integration.
      enum class Dc1cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the DC1C model; deferred until port integration.
      enum class Dc1cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the DC1C model; deferred until port integration.
      enum class Dc1cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the DC1C model; deferred until implementation.
      enum class Dc1cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for DC1C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Dc1cData : public ComponentData<real_type,
                                             index_type,
                                             Dc1cParameters,
                                             Dc1cBuses,
                                             Dc1cSignalInputs,
                                             Dc1cSignalOutputs,
                                             Dc1cMonitorableVariables>
      {
        Dc1cData() = default;

        using Parameters           = Dc1cParameters;
        using Buses                = Dc1cBuses;
        using SignalInputs         = Dc1cSignalInputs;
        using SignalOutputs        = Dc1cSignalOutputs;
        using MonitorableVariables = Dc1cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
