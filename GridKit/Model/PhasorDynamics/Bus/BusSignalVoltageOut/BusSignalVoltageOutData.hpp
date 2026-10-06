/**
 * @file BusSignalVoltageOutData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Signal port definitions for BusSignalVoltageOut.
 *
 * BusSignalVoltageOut is constructed from BusData like other buses. This header adds
 * the port enumerations and a ComponentData alias that satisfies the
 * ModelData concept, so the generic SignalPorts container can be reused for
 * the bus's signal ports and connected from component-style data.
 */
#pragma once

#include <cstddef>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// BusSignalVoltageOut has no bus terminals of its own; it is a bus.
    enum class BusSignalVoltageOutBuses : size_t
    {
    };

    /// Signal inlets of a bus with voltage signal outlets and current signal inlets (see BusSignalVoltageOut)
    enum class BusSignalVoltageOutInputs : size_t
    {
      ir, ///< Real current injection, added to the real current residual
      ii, ///< Imaginary current injection, added to the imaginary current residual
    };

    /// Signal outlets of a bus with voltage signal outlets and current signal inlets (see BusSignalVoltageOut)
    enum class BusSignalVoltageOutOutputs : size_t
    {
      vr, ///< Bus voltage, real component
      vi, ///< Bus voltage, imaginary component
    };

    /**
     * @brief Component-style data for BusSignalVoltageOut signal ports
     *
     * Reuses the bus parameter and monitorable-variable enumerations from
     * BusData. Only the signal maps are used by the bus.
     */
    template <typename real_type, typename index_type>
    using BusSignalVoltageOutData =
        ComponentData<real_type,
                      index_type,
                      BusParameters,
                      BusSignalVoltageOutBuses,
                      BusSignalVoltageOutInputs,
                      BusSignalVoltageOutOutputs,
                      BusMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
