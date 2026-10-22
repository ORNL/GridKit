/**
 * @file DefiniteTimeRelayData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the DefiniteTimeRelay model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Parameter keys for a DefiniteTimeRelay. Mpickup, Mref, and Ttrip are
      /// required; Lockout retains its documented default when omitted.
      enum class DefiniteTimeRelayParameters : size_t
      {
        Mpickup, ///< \f$M_{\mathrm{pickup}}\f$ Pickup level [p.u.]
        Mref,    ///< \f$M_{\mathrm{ref}}\f$ Reference level [p.u.]
        Ttrip,   ///< \f$T_{\mathrm{trip}}\f$ Trip time [sec]
        Lockout, ///< \f$s_{\mathrm{lock}}\f$ Lockout selector [boolean]
      };

      /// Buses for a DefiniteTimeRelay.
      enum class DefiniteTimeRelayBuses : size_t
      {
      };

      /// Signal inputs for a DefiniteTimeRelay.
      enum class DefiniteTimeRelaySignalInputs : size_t
      {
        mr, ///< \f$m_r\f$ Required Known measured quantity, real component [p.u.]
        mi, ///< \f$m_i\f$ Optional Known measured quantity, imaginary component [p.u.]
      };

      /// Signal outputs for a DefiniteTimeRelay.
      enum class DefiniteTimeRelaySignalOutputs : size_t
      {
        trip, ///< \f$s\f$ Required Known trip command [-]
      };

      /// Variables available through the monitor interface.
      enum class DefiniteTimeRelayMonitorableVariables : size_t
      {
        m,    ///< \f$m\f$ Measured quantity magnitude [p.u.]
        x,    ///< \f$x\f$ Relay latch [-]
        trip, ///< \f$s\f$ Trip command [-]
      };

      /**
       * @brief Model data for a DefiniteTimeRelay: parameter values, the
       *        measured-quantity inputs, the trip output, and
       *        monitored-variable selections.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see DefiniteTimeRelay
       */
      template <typename real_type, typename index_type>
      using DefiniteTimeRelayData =
          ComponentData<real_type,
                        index_type,
                        DefiniteTimeRelayParameters,
                        DefiniteTimeRelayBuses,
                        DefiniteTimeRelaySignalInputs,
                        DefiniteTimeRelaySignalOutputs,
                        DefiniteTimeRelayMonitorableVariables>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
