/**
 * @file SexsPtiData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the SEXS-PTI exciter.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /// Parameter keys for the SEXS-PTI exciter model.
      enum class SexsPtiParameters
      {
        Ta,     ///< Numerator time constant of lag-lead block
        Tb,     ///< Denominator time constant of lag-lead block
        Te,     ///< Exciter field time constant
        K,      ///< Voltage regulator gain
        Efdmax, ///< Maximum excitation output
        Efdmin  ///< Minimum excitation output
      };

      /// Buses for the SEXS-PTI exciter model.
      enum class SexsPtiBuses : size_t
      {
        bus, ///< Unique ID of the terminal bus
        SIZE
      };

      /// Signal inputs for the SEXS-PTI exciter model.
      enum class SexsPtiSignalInputs : size_t
      {
        vref, ///< Unique ID of the optional voltage reference signal
        vs,   ///< Unique ID of the optional stabilizer output signal
        vuel, ///< Unique ID of the optional under-excitation limiter signal
        voel, ///< Unique ID of the optional over-excitation limiter signal
        SIZE
      };

      /// Signal outputs for the SEXS-PTI exciter model.
      enum class SexsPtiSignalOutputs : size_t
      {
        efd, ///< Unique ID of the output efd signal
        SIZE
      };

      /// Monitorable variables for the SEXS-PTI exciter model.
      enum class SexsPtiMonitorableVariables
      {
        efd ///< Field voltage output
      };

      template <typename real_type, typename index_type>
      struct SexsPtiData : public ComponentData<real_type,
                                                index_type,
                                                SexsPtiParameters,
                                                SexsPtiBuses,
                                                SexsPtiSignalInputs,
                                                SexsPtiSignalOutputs,
                                                SexsPtiMonitorableVariables>
      {
        SexsPtiData() = default;

        using Parameters           = SexsPtiParameters;
        using Buses                = SexsPtiBuses;
        using SignalInputs         = SexsPtiSignalInputs;
        using SignalOutputs        = SexsPtiSignalOutputs;
        using MonitorableVariables = SexsPtiMonitorableVariables;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
