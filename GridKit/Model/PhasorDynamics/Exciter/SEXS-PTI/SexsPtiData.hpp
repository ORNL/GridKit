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

      /// Terminals for the SEXS-PTI exciter model.
      enum class SexsPtiTerminals : size_t
      {
        bus, ///< Unique ID of the terminal bus
        SIZE
      };

      /// Input ports for the SEXS-PTI exciter model.
      enum class SexsPtiInputPorts : size_t
      {
        vs, ///< Unique ID of the optional stabilizer output signal
        SIZE
      };

      /// Output ports for the SEXS-PTI exciter model.
      enum class SexsPtiOutputPorts : size_t
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
                                                SexsPtiTerminals,
                                                SexsPtiInputPorts,
                                                SexsPtiOutputPorts,
                                                SexsPtiMonitorableVariables>
      {
        SexsPtiData() = default;

        using Parameters           = SexsPtiParameters;
        using Terminals            = SexsPtiTerminals;
        using InputPorts           = SexsPtiInputPorts;
        using OutputPorts          = SexsPtiOutputPorts;
        using MonitorableVariables = SexsPtiMonitorableVariables;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
