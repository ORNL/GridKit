/**
 * @file SwitchData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT switches
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a switch
    enum class SwitchParameters
    {
      N,    ///< Number of phases
      open, ///< Ganged switch command, true is open and false is closed
    };

    /// Buses for a switch
    enum class SwitchBuses : size_t
    {
      bus1, ///< Unique ID of the terminal 1 bus
      bus2, ///< Unique ID of the terminal 2 bus
      SIZE
    };

    /// Signal inputs supported for a switch
    enum class SwitchSignalInputs : size_t
    {
      SIZE
    };

    /// Signal outputs supported for a switch
    enum class SwitchSignalOutputs : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a switch
    enum class SwitchMonitorableVariables
    {
      open,
      i12a,
      i12b,
      i12c
    };

    /**
     * @brief Contains modeling data for a switch
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct SwitchData : public ComponentData<real_type,
                                             index_type,
                                             SwitchParameters,
                                             SwitchBuses,
                                             SwitchSignalInputs,
                                             SwitchSignalOutputs,
                                             SwitchMonitorableVariables>
    {
      SwitchData() = default;

      using Parameters           = SwitchParameters;
      using Buses                = SwitchBuses;
      using SignalInputs         = SwitchSignalInputs;
      using SignalOutputs        = SwitchSignalOutputs;
      using MonitorableVariables = SwitchMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
