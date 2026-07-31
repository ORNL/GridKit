/**
 * @file Oel4cData.hpp
 * @brief Modeling data for the OEL4C model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the OEL4C model.
      enum class Oel4cParameters
      {
        QREF,   ///< \f$Q_{\mathrm{REF}}\f$ OEL timed reactive power limiter pick up level
        Tdelay, ///< \f$T_{\mathrm{delay}}\f$ OEL integral time constant
        KP,     ///< \f$K_{\mathrm{P}}\f$ OEL proportional gain
        KI,     ///< \f$K_{\mathrm{I}}\f$ OEL integral gain
        Vmin    ///< \f$V_{\mathrm{min}}\f$ OEL minimum output
      };

      /// Bus keys for the OEL4C model; deferred until port integration.
      enum class Oel4cBuses : size_t
      {
        SIZE
      };

      /// Signal input keys for the OEL4C model; deferred until port integration.
      enum class Oel4cSignalInputs : size_t
      {
        SIZE
      };

      /// Signal output keys for the OEL4C model; deferred until port integration.
      enum class Oel4cSignalOutputs : size_t
      {
        SIZE
      };

      /// Monitorable variables for the OEL4C model; deferred until implementation.
      enum class Oel4cMonitorableVariables
      {
        NONE
      };

      /**
       * @brief Model data for OEL4C.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Oel4cData : public ComponentData<real_type,
                                              index_type,
                                              Oel4cParameters,
                                              Oel4cBuses,
                                              Oel4cSignalInputs,
                                              Oel4cSignalOutputs,
                                              Oel4cMonitorableVariables>
      {
        Oel4cData() = default;

        using Parameters           = Oel4cParameters;
        using Buses                = Oel4cBuses;
        using SignalInputs         = Oel4cSignalInputs;
        using SignalOutputs        = Oel4cSignalOutputs;
        using MonitorableVariables = Oel4cMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
