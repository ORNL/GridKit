/**
 * @file BusFaultData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for short-to-ground fault
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Initial parameters for a bus fault
    enum class BusFaultParameters
    {
      state0, ///< Whether or not the fault has happened
      R,      ///< Short to ground resistance
      X,      ///< Short to ground reactance
    };

    /// Terminals supported for a bus fault
    enum class BusFaultTerminals : size_t
    {
      bus, ///< Unique ID of the bus where the fault occurs
      SIZE
    };

    /// Input ports supported for a bus fault
    enum class BusFaultInputPorts : size_t
    {
      control_signal, ///< Unique ID of the bus providing a control signal
      SIZE
    };

    /// Output ports supported for a bus fault
    enum class BusFaultOutputPorts : size_t
    {
      SIZE
    };

    /// Variables able to be monitored for a bus fault
    enum class BusFaultMonitorableVariables
    {
      state,
      ir,
      ii
    };

    /**
     * @brief Contains modeling data for a short-to-ground fault
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct BusFaultData : public ComponentData<real_type,
                                               index_type,
                                               BusFaultParameters,
                                               BusFaultTerminals,
                                               BusFaultInputPorts,
                                               BusFaultOutputPorts,
                                               BusFaultMonitorableVariables>
    {
      BusFaultData() = default;

      using Parameters           = BusFaultParameters;
      using Terminals            = BusFaultTerminals;
      using InputPorts           = BusFaultInputPorts;
      using OutputPorts          = BusFaultOutputPorts;
      using MonitorableVariables = BusFaultMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
