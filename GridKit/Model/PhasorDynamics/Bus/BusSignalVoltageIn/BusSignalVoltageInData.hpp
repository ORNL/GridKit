/**
 * @file BusSignalVoltageInData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Signal port definitions for BusSignalVoltageIn.
 *
 * BusSignalVoltageIn reuses BusData for its parameters, initial values and
 * monitored variables; this header only adds the port enumerations.
 */
#pragma once

#include <cstddef>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Signal input ports of a bus with voltage set by signals (see BusSignalVoltageIn)
    enum class BusSignalVoltageInInputs : size_t
    {
      vr, ///< Bus voltage, real component
      vi, ///< Bus voltage, imaginary component
    };

    /// Signal output ports of a bus with voltage set by signals (see BusSignalVoltageIn)
    enum class BusSignalVoltageInOutputs : size_t
    {
      ir, ///< Sum of real current injections from attached components
      ii, ///< Sum of imaginary current injections from attached components
    };
  } // namespace PhasorDynamics
} // namespace GridKit
