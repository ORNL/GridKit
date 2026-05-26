/**
 * @file IeeestData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for IEEEST Stabilizer
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /**
       * @brief Parameter keys for IEEEST Stabilizer model.
       */
      enum class IeeestParameters
      {
        A1,     ///< Notch filter denominator coefficient
        A2,     ///< Notch filter denominator coefficient
        A3,     ///< Notch filter denominator coefficient
        A4,     ///< Notch filter denominator coefficient
        A5,     ///< Notch filter numerator coefficient
        A6,     ///< Notch filter numerator coefficient
        T1,     ///< Lead-lag 1 numerator time constant
        T2,     ///< Lead-lag 1 denominator time constant
        T3,     ///< Lead-lag 2 numerator time constant
        T4,     ///< Lead-lag 2 denominator time constant
        T5,     ///< Washout numerator time constant
        T6,     ///< Washout denominator time constant
        Ks,     ///< Stabilizer gain
        Lsmin,  ///< Minimum stabilizer output limit
        Lsmax,  ///< Maximum stabilizer output limit
        Vcl,    ///< Lower input cutout threshold (not modeled)
        Vcu,    ///< Upper input cutout threshold (not modeled)
        Tdelay, ///< Input time delay (not modeled)
      };

      /**
       * @brief Port keys for IEEEST Stabilizer model.
       */
      enum class IeeestPorts
      {
        input,  ///< Unique ID of the stabilizer input signal
        output, ///< Unique ID of the stabilizer output signal
      };

      /**
       * @brief Monitorable variables for IEEEST Stabilizer model.
       */
      enum class IeeestMonitorableVariables
      {
        vss, ///< Stabilizer output (limited signal)
      };

      /**
       * @brief Contains modeling data for a IEEEST Stabilizer model.
       *
       * @tparam RealT Real parameter data type
       * @tparam IdxT  Integer parameter data type
       */
      template <typename RealT, typename IdxT>
      struct IeeestData : public ComponentData<RealT,
                                               IdxT,
                                               IeeestParameters,
                                               IeeestPorts,
                                               IeeestMonitorableVariables>
      {
        IeeestData() = default;

        using Parameters           = IeeestParameters;
        using Ports                = IeeestPorts;
        using MonitorableVariables = IeeestMonitorableVariables;
      };

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
