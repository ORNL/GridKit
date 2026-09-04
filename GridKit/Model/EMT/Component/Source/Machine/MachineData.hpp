/**
 * @file MachineData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT synchronous machines
 *
 */
#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Initial parameters for a synchronous machine
    enum class MachineParameters
    {
      N,    ///< Number of phases
      S,    ///< Rated three-phase apparent power
      V,    ///< Rated line-to-line RMS voltage
      f,    ///< Rated electrical frequency
      H,    ///< Inertia constant
      F,    ///< Friction torque factor
      Rs,   ///< Stator winding resistance
      Ll,   ///< Stator leakage inductance
      Lmd,  ///< Unsaturated d-axis magnetizing inductance
      Lmq,  ///< Unsaturated q-axis magnetizing inductance
      L0,   ///< Zero-sequence inductance
      Rfd,  ///< Field winding resistance
      Llfd, ///< Field leakage inductance
      R1d,  ///< d-axis damper resistance
      Ll1d, ///< d-axis damper leakage inductance
      R1q,  ///< q-axis damper 1 resistance
      Ll1q, ///< q-axis damper 1 leakage inductance
      R2q,  ///< q-axis damper 2 resistance
      Ll2q, ///< q-axis damper 2 leakage inductance
      S10,  ///< Saturation factor at 1.0 per unit flux
      S12,  ///< Saturation factor at 1.2 per unit flux
      p0,   ///< Initial active power injection
      q0,   ///< Initial reactive power injection
    };

    /// Inputs supported by a synchronous machine
    enum class MachineInputs : size_t
    {
      bus, ///< Component ID of the connected bus
      pm,  ///< Mechanical-power signal ID from a governor
      efd, ///< Field-voltage signal ID from an exciter
      SIZE
    };

    /// Outputs supported by a synchronous machine
    enum class MachineOutputs : size_t
    {
      speed, ///< Rotor-speed output signal ID
      SIZE
    };

    /// Variables able to be monitored for a synchronous machine
    enum class MachineMonitorableVariables
    {
      theta,
      omega,
      te,
      ifd,
      efd,
      ks,
      psi_at,
      ia,
      ib,
      ic,
      p,
      q,
      id,
      iq
    };

    /**
     * @brief Contains modeling data for a synchronous machine
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct MachineData : public ComponentData<real_type,
                                              index_type,
                                              MachineParameters,
                                              MachineInputs,
                                              MachineOutputs,
                                              MachineMonitorableVariables>
    {
      MachineData() = default;

      using Parameters           = MachineParameters;
      using Inputs               = MachineInputs;
      using Outputs              = MachineOutputs;
      using MonitorableVariables = MachineMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
