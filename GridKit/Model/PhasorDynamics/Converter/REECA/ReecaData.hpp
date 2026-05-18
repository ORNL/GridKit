/**
 * @file ReecaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REECA converter-control model.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Parameter keys for the REECA converter-control model.
      enum class ReecaParameters
      {
        Sbase, ///< REECA model MVA base
        spf,   ///< Power-factor control flag
        sV,    ///< Voltage-control flag
        sQ,    ///< Reactive-power-control flag
        sP,    ///< Active-power-control speed-feedback flag
        sPQ,   ///< P/Q priority flag
        Trv,   ///< Terminal-voltage transducer time constant
        Tp,    ///< Active-power transducer time constant
        Vref0, ///< Optional initial voltage reference
        Vdip,  ///< Low-voltage active-current management threshold
        Vup,   ///< Over-voltage active-current management threshold
        Dbd1,  ///< Lower reactive-current deadband
        Dbd2,  ///< Upper reactive-current deadband
        Kqv,   ///< Voltage-dip reactive-current gain
        Iql1,  ///< Minimum reactive current during voltage dip
        Iqh1,  ///< Maximum reactive current during voltage dip
        Iqfrz, ///< Frozen reactive current during voltage recovery
        Thld,  ///< Low-voltage hold time
        Qmax,  ///< Maximum reactive-power command
        Qmin,  ///< Minimum reactive-power command
        Kqp,   ///< Reactive-power proportional gain
        Kqi,   ///< Reactive-power integral gain
        Vmax,  ///< Maximum voltage-control command
        Vmin,  ///< Minimum voltage-control command
        Vref1, ///< Fixed voltage-reference bias
        Kvp,   ///< Voltage-control proportional gain
        Kvi,   ///< Voltage-control integral gain
        Tiq,   ///< Reactive-current lag time constant
        Tpord, ///< Active-power order lag time constant
        dPmax, ///< Maximum active-power order ramp rate
        dPmin, ///< Minimum active-power order ramp rate
        Pmax,  ///< Maximum active-power command
        Pmin,  ///< Minimum active-power command
        Imax,  ///< Converter current magnitude limit
        vq1,   ///< Reactive-current VDL voltage breakpoint 1
        lq1,   ///< Reactive-current VDL limit 1
        vq2,   ///< Reactive-current VDL voltage breakpoint 2
        lq2,   ///< Reactive-current VDL limit 2
        vq3,   ///< Reactive-current VDL voltage breakpoint 3
        lq3,   ///< Reactive-current VDL limit 3
        vq4,   ///< Reactive-current VDL voltage breakpoint 4
        lq4,   ///< Reactive-current VDL limit 4
        vp1,   ///< Active-current VDL voltage breakpoint 1
        lp1,   ///< Active-current VDL limit 1
        vp2,   ///< Active-current VDL voltage breakpoint 2
        lp2,   ///< Active-current VDL limit 2
        vp3,   ///< Active-current VDL voltage breakpoint 3
        lp3,   ///< Active-current VDL limit 3
        vp4,   ///< Active-current VDL voltage breakpoint 4
        lp4,   ///< Active-current VDL limit 4
        Thld2  ///< Post-recovery active-current hold time
      };

      /// Ports for the REECA converter-control model.
      enum class ReecaPorts
      {
        bus,    ///< Terminal bus ID
        pe,     ///< Required active-power feedback signal ID
        qgen,   ///< Required reactive-power feedback signal ID
        omega,  ///< Optional speed-deviation signal ID
        qext,   ///< Optional external reactive-power command signal ID
        pfaref, ///< Power-factor angle reference signal ID when spf is true
        pref,   ///< Optional external active-power reference signal ID
        iqcmd,  ///< Optional reactive-current command output signal ID
        ipcmd   ///< Optional active-current command output signal ID
      };

      /// External variables of a `Reeca`.
      enum class ReecaExternalVariables : std::size_t
      {
        PE,     ///< Required active-power feedback signal
        QGEN,   ///< Required reactive-power feedback signal
        OMEGA,  ///< Optional speed-deviation signal
        QEXT,   ///< Optional external reactive-power command signal
        PFAREF, ///< Power-factor angle reference signal
        PREF,   ///< Optional active-power reference signal
        MAXIMUM,
      };

      /// Variables available through the monitor interface.
      enum class ReecaMonitorableVariables
      {
        iqcmd,  ///< Reactive-current command output
        ipcmd,  ///< Active-current command output
        vmeas,  ///< Filtered terminal-voltage measurement
        pmeas,  ///< Filtered active-power measurement
        piq,    ///< Reactive-power PI controller state
        piv,    ///< Voltage-control PI controller state
        qv,     ///< Reactive-current lag state
        pord,   ///< Active-power order state
        qref,   ///< Reactive-power reference placeholder
        sdip,   ///< Voltage dip/overvoltage indicator placeholder
        iqmax,  ///< Reactive-current maximum placeholder
        ipmax,  ///< Active-current maximum placeholder
        iqv,    ///< Voltage-support reactive-current command placeholder
        vqctrl, ///< Reactive-power-control output placeholder
        iqbase  ///< Base reactive-current command placeholder
      };

      template <typename RealT, typename IdxT>
      struct ReecaData : public ComponentData<RealT,
                                              IdxT,
                                              ReecaParameters,
                                              ReecaPorts,
                                              ReecaMonitorableVariables,
                                              ReecaExternalVariables>
      {
        ReecaData() = default;

        using Parameters           = ReecaParameters;
        using Ports                = ReecaPorts;
        using MonitorableVariables = ReecaMonitorableVariables;
        using ExternalVariables    = ReecaExternalVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
