/**
 * @file Oel3cData.hpp
 * @brief Modeling data for the OEL3C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the OEL3C model.
      enum class Oel3cParameters
      {
        ITFpu,    ///< \f$I_{\mathrm{TFpu}}\f$ OEL timed field current limiter pick up level
        KSCALE,   ///< \f$K_{\mathrm{SCALE}}\f$ OEL input signal scaling factor
        TF,       ///< \f$T_{\mathrm{F}}\f$ OEL field current measurement time constant
        K1,       ///< \f$K_{1}\f$ Exponent for OEL error calculation
        KOEL,     ///< \f$K_{\mathrm{OEL}}\f$ OEL gain
        TOEL,     ///< \f$T_{\mathrm{OEL}}\f$ OEL integral time constant
        KPOEL,    ///< \f$K_{\mathrm{POEL}}\f$ OEL proportional gain
        VOELmax1, ///< \f$V_{\mathrm{OELmax1}}\f$ OEL integrator maximum output
        VOELmin1, ///< \f$V_{\mathrm{OELmin1}}\f$ OEL integrator minimum output
        VOELmax2, ///< \f$V_{\mathrm{OELmax2}}\f$ OEL maximum output
        VOELmin2  ///< \f$V_{\mathrm{OELmin2}}\f$ OEL minimum output
      };

      /// Bus keys for the OEL3C model; deferred until port integration.
      enum class Oel3cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the OEL3C model; deferred until port integration.
      enum class Oel3cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the OEL3C model; deferred until port integration.
      enum class Oel3cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the OEL3C model; deferred until implementation.
      enum class Oel3cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for OEL3C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Oel3cData : public ComponentData<real_type,
                                              index_type,
                                              Oel3cParameters,
                                              Oel3cBuses,
                                              Oel3cSignalInputs,
                                              Oel3cSignalOutputs,
                                              Oel3cMonitorableVariables>
      {
        Oel3cData() = default;

        using Parameters           = Oel3cParameters;
        using Buses                = Oel3cBuses;
        using SignalInputs         = Oel3cSignalInputs;
        using SignalOutputs        = Oel3cSignalOutputs;
        using MonitorableVariables = Oel3cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
