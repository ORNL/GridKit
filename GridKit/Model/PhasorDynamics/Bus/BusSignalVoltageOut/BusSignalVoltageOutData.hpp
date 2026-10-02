/**
 * @file BusSignalVoltageOutData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Signal port definitions for BusSignalVoltageOut.
 *
 * BusSignalVoltageOut reuses BusData for its parameters, initial values and
 * monitored variables; this header only adds the port enumerations.
 */
#pragma once

#include <cstddef>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Signal input ports of a bus with signal ports (see BusSignalVoltageOut)
    enum class BusSignalVoltageOutInputs : size_t
    {
      ir, ///< Real current injection, added to the real current residual
      ii, ///< Imaginary current injection, added to the imaginary current residual
    };

    /// Signal output ports of a bus with signal ports (see BusSignalVoltageOut)
    enum class BusSignalVoltageOutOutputs : size_t
    {
      vr, ///< Bus voltage, real component
      vi, ///< Bus voltage, imaginary component
    };
  } // namespace PhasorDynamics
} // namespace GridKit
