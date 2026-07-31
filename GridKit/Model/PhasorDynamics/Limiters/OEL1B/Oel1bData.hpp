/**
 * @file Oel1bData.hpp
 * @brief Modeling data for the OEL1B overexcitation-limiter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Limiter
    {
      /// Parameter keys for the OEL1B overexcitation-limiter model.
      enum class Oel1bParameters
      {
        Itfpu,   ///< \f$I_{\mathrm{TFPU}}\f$ Timed limiter pickup level
        Ifdmax,  ///< \f$I_{\mathrm{FDMAX}}\f$ Instantaneous field-current limit
        Ifdlim,  ///< \f$I_{\mathrm{FDLIM}}\f$ Timed field-current limit
        Hyst,    ///< \f$\mathrm{HYST}\f$ Pickup and dropout hysteresis
        Kcd,     ///< \f$K_{\mathrm{CD}}\f$ Cool-down gain
        Kramp,   ///< \f$K_{\mathrm{RAMP}}\f$ Limiter output ramp rate
        Ifdrated ///< \f$I_{\mathrm{FDrated}}\f$ Rated field current
      };

      /// Buses for the OEL1B overexcitation-limiter model.
      enum class Oel1bBuses : size_t
      {
        SIZE
      };

      /// Signal inputs for the OEL1B overexcitation-limiter model.
      enum class Oel1bSignalInputs : size_t
      {
        ifd, ///< Required generator field-current signal ID
        SIZE
      };

      /// Signal outputs for the OEL1B overexcitation-limiter model.
      enum class Oel1bSignalOutputs : size_t
      {
        voel, ///< Required overexcitation-limiter output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class Oel1bMonitorableVariables
      {
        voel,  ///< Overexcitation-limiter output
        timer, ///< Accumulated overexcitation time
        hyst   ///< Active hysteresis value
      };

      /**
       * @brief Model data for OEL1B: parameters, signal connections, and
       *        monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       */
      template <typename real_type, typename index_type>
      struct Oel1bData : public ComponentData<real_type,
                                              index_type,
                                              Oel1bParameters,
                                              Oel1bBuses,
                                              Oel1bSignalInputs,
                                              Oel1bSignalOutputs,
                                              Oel1bMonitorableVariables>
      {
        Oel1bData() = default;

        using Parameters           = Oel1bParameters;
        using Buses                = Oel1bBuses;
        using SignalInputs         = Oel1bSignalInputs;
        using SignalOutputs        = Oel1bSignalOutputs;
        using MonitorableVariables = Oel1bMonitorableVariables;
      };
    } // namespace Limiter
  } // namespace PhasorDynamics
} // namespace GridKit
