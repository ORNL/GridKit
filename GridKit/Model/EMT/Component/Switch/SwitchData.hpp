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

    /// Inputs supported by a switch
    enum class SwitchInputs : size_t
    {
      bus1, ///< Component ID of the terminal 1 bus
      bus2, ///< Component ID of the terminal 2 bus
      SIZE
    };

    /// Outputs supported by a switch
    enum class SwitchOutputs : size_t
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
                                             SwitchInputs,
                                             SwitchOutputs,
                                             SwitchMonitorableVariables>
    {
      SwitchData() = default;

      using Parameters           = SwitchParameters;
      using Inputs               = SwitchInputs;
      using Outputs              = SwitchOutputs;
      using MonitorableVariables = SwitchMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
