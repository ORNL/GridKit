/**
 * @file ForcedOscillationData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the forced oscillation signal operator.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalOperator
    {
      /// Parameter keys for the forced oscillation signal operator.
      enum class ForcedOscillationParameters
      {
        A,    ///< Oscillation amplitude
        f,    ///< Oscillation frequency
        Kf,   ///< Linear frequency ramp
        Phi,  ///< Phase offset
        Bias, ///< DC offset
        Kin,  ///< Input gain
        u0,   ///< Standalone input value
        Ton,  ///< Start time
        Toff, ///< Stop time
        Tr,   ///< Rise time
        Tf,   ///< Fall time
        Kd,   ///< Exponential decay rate
        Lmin, ///< Output lower limit
        Lmax  ///< Output upper limit
      };

      /// Port keys for the forced oscillation signal operator.
      enum class ForcedOscillationPorts
      {
        input, ///< Optional input signal ID
        output ///< Output signal ID
      };

      /// Variables available through the monitor interface.
      enum class ForcedOscillationMonitorableVariables
      {
        in,    ///< Input value
        env,   ///< Envelope value
        force, ///< Oscillatory forcing
        vraw,  ///< Unlimited output
        out,   ///< Limited output
        active ///< Active-window indicator
      };

      /**
       * @brief Contains modeling data for a forced oscillation signal operator.
       *
       * @tparam RealT Real parameter data type
       * @tparam IdxT  Integer parameter data type
       */
      template <typename RealT, typename IdxT>
      struct ForcedOscillationData : public ComponentData<RealT,
                                                          IdxT,
                                                          ForcedOscillationParameters,
                                                          ForcedOscillationPorts,
                                                          ForcedOscillationMonitorableVariables>
      {
        ForcedOscillationData() = default;

        using Parameters           = ForcedOscillationParameters;
        using Ports                = ForcedOscillationPorts;
        using MonitorableVariables = ForcedOscillationMonitorableVariables;
      };

    } // namespace SignalOperator
  } // namespace PhasorDynamics
} // namespace GridKit
