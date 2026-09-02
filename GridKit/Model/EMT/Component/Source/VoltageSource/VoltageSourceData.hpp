/**
 * @file VoltageSourceData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT voltage sources
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a voltage source
    enum class VoltageSourceParameters
    {
      N,     ///< Number of phases
      E,     ///< Source voltage magnitudes, RMS
      phi,   ///< Source phase offsets
      omega, ///< Source angular frequency
      Rs,    ///< Series resistance matrix
      Ls,    ///< Series inductance matrix
    };

    /// Buses for a voltage source
    enum class VoltageSourceBuses : size_t
    {
      bus, ///< Unique ID of the bus to which the source is connected
      SIZE
    };

    /// Signal inputs supported for a voltage source
    enum class VoltageSourceSignalInputs : size_t
    {
      SIZE
    };

    /// Signal outputs supported for a voltage source
    enum class VoltageSourceSignalOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a voltage source
    enum class VoltageSourceMonitorableVariables
    {
      ea,
      eb,
      ec,
      ia,
      ib,
      ic
    };

    /**
     * @brief Contains modeling data for a voltage source
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct VoltageSourceData : public ComponentData<real_type,
                                                    index_type,
                                                    VoltageSourceParameters,
                                                    VoltageSourceBuses,
                                                    VoltageSourceSignalInputs,
                                                    VoltageSourceSignalOutputs,
                                                    VoltageSourceMonitorableVariables>
    {
      VoltageSourceData() = default;

      using Parameters           = VoltageSourceParameters;
      using Buses                = VoltageSourceBuses;
      using SignalInputs         = VoltageSourceSignalInputs;
      using SignalOutputs        = VoltageSourceSignalOutputs;
      using MonitorableVariables = VoltageSourceMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
