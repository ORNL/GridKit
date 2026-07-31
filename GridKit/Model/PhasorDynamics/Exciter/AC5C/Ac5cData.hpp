/**
 * @file Ac5cData.hpp
 * @brief Modeling data for the AC5C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC5C model.
      enum class Ac5cParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        KA,    ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        TA,    ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        VAmax, ///< \f$V_{\mathrm{Amax}}\f$ Maximum regulator output
        VAmin, ///< \f$V_{\mathrm{Amin}}\f$ Minimum regulator output
        TE,    ///< \f$T_{\mathrm{E}}\f$ Exciter field time constant
        KE,    ///< \f$K_{\mathrm{E}}\f$ Exciter field proportional constant
        E1,    ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        Se1,   ///< \f$S_{\mathrm{E}}(E_{1})\f$ Exciter saturation factor at exciter output voltage E 1
        E2,    ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        Se2,   ///< \f$S_{\mathrm{E}}(E_{2})\f$ Exciter saturation factor at exciter output voltage E 2
        KF,    ///< \f$K_{\mathrm{F}}\f$ Rate feedback excitation system stabilizer gain
        TF1,   ///< \f$T_{\mathrm{F1}}\f$ Rate feedback excitation system stabilizer time constant
        TF2,   ///< \f$T_{\mathrm{F2}}\f$ Rate feedback excitation system stabilizer time constant
        TF3    ///< \f$T_{\mathrm{F3}}\f$ Rate feedback excitation system stabilizer time constant
      };

      /// Bus keys for the AC5C model; deferred until port integration.
      enum class Ac5cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC5C model; deferred until port integration.
      enum class Ac5cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC5C model; deferred until port integration.
      enum class Ac5cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC5C model; deferred until implementation.
      enum class Ac5cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC5C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac5cData : public ComponentData<real_type,
                                             index_type,
                                             Ac5cParameters,
                                             Ac5cBuses,
                                             Ac5cSignalInputs,
                                             Ac5cSignalOutputs,
                                             Ac5cMonitorableVariables>
      {
        Ac5cData() = default;

        using Parameters           = Ac5cParameters;
        using Buses                = Ac5cBuses;
        using SignalInputs         = Ac5cSignalInputs;
        using SignalOutputs        = Ac5cSignalOutputs;
        using MonitorableVariables = Ac5cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
