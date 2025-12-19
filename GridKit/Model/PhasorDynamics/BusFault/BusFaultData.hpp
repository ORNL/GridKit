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

    /// Ports supported for a bus fault
    enum class BusFaultPorts
    {
      bus,            ///< Unique ID of the bus where the fault occurs
      control_signal, ///< Unique ID of the bus providing a control signal
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
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename RealT, typename IdxT>
    struct BusFaultData : public ComponentData<RealT,
                                               IdxT,
                                               BusFaultParameters,
                                               BusFaultPorts,
                                               BusFaultMonitorableVariables>
    {
      BusFaultData() = default;

      using Parameters           = BusFaultParameters;
      using Ports                = BusFaultPorts;
      using MonitorableVariables = BusFaultMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
