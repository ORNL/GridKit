/**
 * @file Dc2cData.hpp
 * @brief Modeling data for the DC2C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the DC2C model.
      enum class Dc2cParameters
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
        VRmax, ///< \f$V_{\mathrm{Rmax}}\f$ Maximum controller output
        VRmin, ///< \f$V_{\mathrm{Rmin}}\f$ Minimum controller output
        KF,    ///< \f$K_{\mathrm{F}}\f$ Rate feedback gain
        TF,    ///< \f$T_{\mathrm{F}}\f$ Rate feedback time constant
        E1,    ///< \f$E_{1}\f$ Exciter output voltage for saturation factor S (E ) E 1
        SE1,   ///< \f$S_{\mathrm{E1}}\f$ Exciter saturation factor at exciter output voltage E 1
        E2,    ///< \f$E_{2}\f$ Exciter output voltage for saturation factor S (E ) E 2
        SE2    ///< \f$S_{\mathrm{E2}}\f$ Exciter saturation factor at exciter output voltage E 2
      };

      /// Bus keys for the DC2C model; deferred until port integration.
      enum class Dc2cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the DC2C model; deferred until port integration.
      enum class Dc2cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the DC2C model; deferred until port integration.
      enum class Dc2cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the DC2C model; deferred until implementation.
      enum class Dc2cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for DC2C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Dc2cData : public ComponentData<real_type,
                                             index_type,
                                             Dc2cParameters,
                                             Dc2cBuses,
                                             Dc2cSignalInputs,
                                             Dc2cSignalOutputs,
                                             Dc2cMonitorableVariables>
      {
        Dc2cData() = default;

        using Parameters           = Dc2cParameters;
        using Buses                = Dc2cBuses;
        using SignalInputs         = Dc2cSignalInputs;
        using SignalOutputs        = Dc2cSignalOutputs;
        using MonitorableVariables = Dc2cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
