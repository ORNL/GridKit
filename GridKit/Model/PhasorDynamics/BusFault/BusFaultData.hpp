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
      R, ///< \f$R\f$ Short to ground resistance
      X, ///< \f$X\f$ Short to ground reactance
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
      status, ///< Unique ID of the fault status signal
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
      status,
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
                                               BusFaultBuses,
                                               BusFaultSignalInputs,
                                               BusFaultSignalOutputs,
                                               BusFaultMonitorableVariables>
    {
      BusFaultData() = default;

      using Parameters           = BusFaultParameters;
      using Buses                = BusFaultBuses;
      using SignalInputs         = BusFaultSignalInputs;
      using SignalOutputs        = BusFaultSignalOutputs;
      using MonitorableVariables = BusFaultMonitorableVariables;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
