/**
 * @file OvercurrentRelayData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the overcurrent relay model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Parameter keys for the overcurrent relay model. Both parameters are
      /// required relay settings.
      enum class OvercurrentRelayParameters : size_t
      {
        Ipickup, ///< \f$I_{\mathrm{pickup}}\f$ Pickup current magnitude on system base [p.u.]
        Ttrip,   ///< \f$T_{\mathrm{trip}}\f$ Time from sustained pickup to trip [sec]
      };

      /// Buses for the overcurrent relay model.
      enum class OvercurrentRelayBuses : size_t
      {
      };

      /// Signal inputs for the overcurrent relay model.
      enum class OvercurrentRelaySignalInputs : size_t
      {
        ir, ///< \f$I_{\mathrm{r}}\f$ Required Known measured-current real-component input [p.u.]
        ii, ///< \f$I_{\mathrm{i}}\f$ Required Known measured-current imaginary-component input [p.u.]
      };

      /// Signal outputs for the overcurrent relay model.
      enum class OvercurrentRelaySignalOutputs : size_t
      {
        trip, ///< \f$s\f$ Required Known trip-command output [-]
      };

      /// Variables available through the monitor interface.
      enum class OvercurrentRelayMonitorableVariables : size_t
      {
        im,   ///< \f$|I|\f$ Measured current magnitude [p.u.]
        x,    ///< \f$x\f$ Lockout latch state [-]
        trip, ///< \f$s\f$ Trip command [-]
      };

      /**
       * @brief Model data for overcurrent relay parameters, signal ports, and
       *        monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see OvercurrentRelay
       */
      template <typename real_type, typename index_type>
      using OvercurrentRelayData =
          ComponentData<real_type,
                        index_type,
                        OvercurrentRelayParameters,
                        OvercurrentRelayBuses,
                        OvercurrentRelaySignalInputs,
                        OvercurrentRelaySignalOutputs,
                        OvercurrentRelayMonitorableVariables>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
