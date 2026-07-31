/**
 * @file Ac4cData.hpp
 * @brief Modeling data for the AC4C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the AC4C model.
      enum class Ac4cParameters
      {
        RC,    ///< \f$R_{\mathrm{C}}\f$ Resistive component of load compensation
        XC,    ///< \f$X_{\mathrm{C}}\f$ Reactance component of load compensation
        TR,    ///< \f$T_{\mathrm{R}}\f$ Regulator input filter time constant
        TB,    ///< \f$T_{\mathrm{B}}\f$ Regulator denominator (lag) time constant
        TC,    ///< \f$T_{\mathrm{C}}\f$ Regulator numerator (lead) time constant
        TA,    ///< \f$T_{\mathrm{A}}\f$ Regulator output time constant
        VImax, ///< \f$V_{\mathrm{Imax}}\f$ Voltage regulator input (voltage error) maximum limit
        VImin, ///< \f$V_{\mathrm{Imin}}\f$ Voltage regulator input (voltage error) minimum limit
        VRmax, ///< \f$V_{\mathrm{Rmax}}\f$ Maximum regulator output
        VRmin, ///< \f$V_{\mathrm{Rmin}}\f$ Minimum regulator output
        KA,    ///< \f$K_{\mathrm{A}}\f$ Regulator output gain
        KC     ///< \f$K_{\mathrm{C}}\f$ Rectifier loading factor proportional to commutating reactance
      };

      /// Bus keys for the AC4C model; deferred until port integration.
      enum class Ac4cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the AC4C model; deferred until port integration.
      enum class Ac4cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the AC4C model; deferred until port integration.
      enum class Ac4cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the AC4C model; deferred until implementation.
      enum class Ac4cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for AC4C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Ac4cData : public ComponentData<real_type,
                                             index_type,
                                             Ac4cParameters,
                                             Ac4cBuses,
                                             Ac4cSignalInputs,
                                             Ac4cSignalOutputs,
                                             Ac4cMonitorableVariables>
      {
        Ac4cData() = default;

        using Parameters           = Ac4cParameters;
        using Buses                = Ac4cBuses;
        using SignalInputs         = Ac4cSignalInputs;
        using SignalOutputs        = Ac4cSignalOutputs;
        using MonitorableVariables = Ac4cMonitorableVariables;
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
