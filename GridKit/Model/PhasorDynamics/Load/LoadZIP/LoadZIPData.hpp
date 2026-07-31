
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Parameters for a loadZIP
    enum class LoadZIPParameters
    {
      Pnom,   ///< Legacy fallback for the p input
      Qnom,   ///< Legacy fallback for the q input
      Vnom,   ///< Nominal voltage magnitude
      alphaI, ///< Fraction of load represented as constant current
      alphaP, ///< Fraction of load represented as constant power
    };

    /// Buses for a loadZIP
    enum class LoadZIPBuses : size_t
    {
      bus, ///< Unique ID of the bus to which the loadZIP is connected
      SIZE
    };

    /// Signal inputs supported for a loadZIP
    enum class LoadZIPSignalInputs : size_t
    {
      p,      ///< Initial terminal active-power injection
      q,      ///< Initial terminal reactive-power injection
      online, ///< Connection status (zero is disconnected, nonzero is connected)
      SIZE
    };

    /// Signal outputs supported for a loadZIP
    enum class LoadZIPSignalOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a loadZIP
    enum class LoadZIPMonitorableVariables
    {
      ir,
      ii,
      im,
      p,
      q
    };

    /**
     * @brief Contains modeling data for a load
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct LoadZIPData : public ComponentData<real_type,
                                              index_type,
                                              LoadZIPParameters,
                                              LoadZIPBuses,
                                              LoadZIPSignalInputs,
                                              LoadZIPSignalOutputs,
                                              LoadZIPMonitorableVariables>
    {
      LoadZIPData() = default;

      using Parameters           = LoadZIPParameters;
      using Buses                = LoadZIPBuses;
      using SignalInputs         = LoadZIPSignalInputs;
      using SignalOutputs        = LoadZIPSignalOutputs;
      using MonitorableVariables = LoadZIPMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
