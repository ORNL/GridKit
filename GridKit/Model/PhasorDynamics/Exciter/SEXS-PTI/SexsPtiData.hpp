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

      /// Ports for the SEXS-PTI exciter model.
      enum class SexsPtiPorts
      {
        bus, ///< Unique ID of the terminal bus
        efd, ///< Unique ID of the output efd signal
        vs   ///< Unique ID of the optional stabilizer output signal
      };

      /// Monitorable variables for the SEXS-PTI exciter model.
      enum class SexsPtiMonitorableVariables
      {
        efd ///< Field voltage output
      };

      template <typename RealT, typename IdxT>
      struct SexsPtiData : public ComponentData<RealT,
                                                IdxT,
                                                SexsPtiParameters,
                                                SexsPtiPorts,
                                                SexsPtiMonitorableVariables>
      {
        SexsPtiData() = default;

        using Parameters           = SexsPtiParameters;
        using Ports                = SexsPtiPorts;
        using MonitorableVariables = SexsPtiMonitorableVariables;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
