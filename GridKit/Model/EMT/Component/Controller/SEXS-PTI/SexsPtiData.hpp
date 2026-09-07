/**
 * @file SexsPtiData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the SEXS-PTI exciter.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Parameter keys for the SEXS-PTI exciter model.
      enum class SexsPtiParameters
      {
        V,      ///< Rated line-to-line RMS terminal voltage in volts
        Tr,     ///< Optional terminal-voltage measurement lag in seconds (zero bypasses)
        Ta,     ///< Numerator time constant of lag-lead block
        Tb,     ///< Denominator time constant of lag-lead block
        Te,     ///< Exciter field time constant
        K,      ///< Voltage regulator gain
        Efdmax, ///< Maximum excitation output
        Efdmin  ///< Minimum excitation output
      };

      /// Signal inputs for the SEXS-PTI exciter model.
      enum class SexsPtiInputs : size_t
      {
        va,   ///< Phase-a terminal voltage
        vb,   ///< Phase-b terminal voltage
        vc,   ///< Phase-c terminal voltage
        vref, ///< Unique ID of the optional voltage reference signal
        vs,   ///< Unique ID of the optional stabilizer output signal
        vuel, ///< Unique ID of the optional under-excitation limiter signal
        voel, ///< Unique ID of the optional over-excitation limiter signal
        SIZE
      };

      /// Signal outputs for the SEXS-PTI exciter model.
      enum class SexsPtiOutputs : size_t
      {
        efd, ///< Unique ID of the output efd signal
        SIZE
      };

      /// Monitorable variables for the SEXS-PTI exciter model.
      enum class SexsPtiMonitorableVariables
      {
        efd, ///< Field voltage output
        vts, ///< Measured terminal voltage
        vr,  ///< Lead-lag state
        vtr  ///< Terminal voltage error
      };

      template <typename real_type, typename index_type>
      struct SexsPtiData : public ComponentData<real_type,
                                                index_type,
                                                SexsPtiParameters,
                                                SexsPtiInputs,
                                                SexsPtiOutputs,
                                                SexsPtiMonitorableVariables>
      {
        SexsPtiData() = default;

        using Parameters           = SexsPtiParameters;
        using Inputs               = SexsPtiInputs;
        using Outputs              = SexsPtiOutputs;
        using MonitorableVariables = SexsPtiMonitorableVariables;
      };

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
