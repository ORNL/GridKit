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

    /// Buses supported for a bus fault
    enum class BusFaultBuses : size_t
    {
      bus, ///< Unique ID of the bus where the fault occurs
      SIZE
    };

    /// Signal inputs supported for a bus fault
    enum class BusFaultSignalInputs : size_t
    {
      control_signal, ///< Unique ID of the bus providing a control signal
      SIZE
    };

    /// Signal outputs supported for a bus fault
    enum class BusFaultSignalOutputs : size_t
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
    using BusFaultData =
        ComponentData<real_type,
                      index_type,
                      BusFaultParameters,
                      BusFaultBuses,
                      BusFaultSignalInputs,
                      BusFaultSignalOutputs,
                      BusFaultMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
