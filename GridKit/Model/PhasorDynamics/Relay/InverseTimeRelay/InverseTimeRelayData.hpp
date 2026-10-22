/**
 * @file InverseTimeRelayData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the InverseTimeRelay model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /// Parameter keys for an InverseTimeRelay. Ipickup, TD, A, B, and alpha
      /// are required; tr and Lockout retain their documented defaults when
      /// omitted.
      enum class InverseTimeRelayParameters : size_t
      {
        Ipickup, ///< \f$I_{\mathrm{pickup}}\f$ Pickup current [p.u.]
        TD,      ///< \f$D\f$ Time dial [-]
        A,       ///< \f$A\f$ Curve coefficient [sec]
        B,       ///< \f$B\f$ Curve offset [sec]
        alpha,   ///< \f$\alpha\f$ Curve exponent [-]
        tr,      ///< \f$t_r\f$ Reset constant [sec]
        Lockout, ///< \f$s_{\mathrm{lock}}\f$ Lockout selector [boolean]
      };

      /// Buses for an InverseTimeRelay.
      enum class InverseTimeRelayBuses : size_t
      {
      };

      /// Signal inputs for an InverseTimeRelay.
      enum class InverseTimeRelaySignalInputs : size_t
      {
        ir, ///< \f$I_r\f$ Required Known measured current, real component [p.u.]
        ii, ///< \f$I_i\f$ Required Known measured current, imaginary component [p.u.]
      };

      /// Signal outputs for an InverseTimeRelay.
      enum class InverseTimeRelaySignalOutputs : size_t
      {
        trip, ///< \f$s\f$ Required Known trip command [-]
      };

      /// Variables available through the monitor interface.
      enum class InverseTimeRelayMonitorableVariables : size_t
      {
        im,   ///< \f$I_m\f$ Measured current magnitude [p.u.]
        x,    ///< \f$x\f$ Relay latch [-]
        trip, ///< \f$s\f$ Trip command [-]
      };

      /**
       * @brief Model data for an InverseTimeRelay: parameter values, the
       *        measured-current inputs, the trip output, and
       *        monitored-variable selections.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see InverseTimeRelay
       */
      template <typename real_type, typename index_type>
      using InverseTimeRelayData =
          ComponentData<real_type,
                        index_type,
                        InverseTimeRelayParameters,
                        InverseTimeRelayBuses,
                        InverseTimeRelaySignalInputs,
                        InverseTimeRelaySignalOutputs,
                        InverseTimeRelayMonitorableVariables>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
