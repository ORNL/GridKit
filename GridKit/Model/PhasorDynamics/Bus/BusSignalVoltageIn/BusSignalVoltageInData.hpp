/**
 * @file BusSignalVoltageInData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Signal port definitions for BusSignalVoltageIn.
 *
 * BusSignalVoltageIn is constructed from BusData like other buses. This header adds
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
    /// BusSignalVoltageIn has no bus terminals of its own; it is a bus.
    enum class BusSignalVoltageInBuses : size_t
    {
    };

    /// Signal inlets of a bus with voltage signal inlets and current signal outlets (see BusSignalVoltageIn)
    enum class BusSignalVoltageInInputs : size_t
    {
      vr, ///< Bus voltage, real component
      vi, ///< Bus voltage, imaginary component
    };

    /// Signal outlets of a bus with voltage signal inlets and current signal outlets (see BusSignalVoltageIn)
    enum class BusSignalVoltageInOutputs : size_t
    {
      ir, ///< Sum of real current injections from attached components
      ii, ///< Sum of imaginary current injections from attached components
    };

    /**
     * @brief Component-style data for BusSignalVoltageIn signal ports
     *
     * Reuses the bus parameter and monitorable-variable enumerations from
     * BusData. Only the signal maps are used by the bus.
     */
    template <typename real_type, typename index_type>
    using BusSignalVoltageInData =
        ComponentData<real_type,
                      index_type,
                      BusParameters,
                      BusSignalVoltageInBuses,
                      BusSignalVoltageInInputs,
                      BusSignalVoltageInOutputs,
                      BusMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
