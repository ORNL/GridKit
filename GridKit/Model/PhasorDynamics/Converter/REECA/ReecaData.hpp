/**
 * @file ReecaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REECA electrical-control model.
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
      /// Parameter keys for the REECA electrical-control model.
      enum class ReecaParameters
      {
        mva,    ///< Model MVA base
        PfFlag, ///< Power-factor control flag
        VFlag,  ///< Voltage-control mode flag
        QFlag,  ///< Reactive-power control flag
        PFlag,  ///< Active-power reference speed-proxy multiplier flag
        Pqflag, ///< P/Q current-priority flag
        Trv,    ///< Voltage-measurement filter time constant
        Tp,     ///< Electrical-power measurement filter time constant
        Vref0,  ///< Outer-loop voltage reference
        Vdip,   ///< Low-voltage dip threshold
        Vup,    ///< High-voltage threshold
        dbd1,   ///< Lower voltage-error deadband threshold
        dbd2,   ///< Upper voltage-error deadband threshold
        kqv,    ///< Reactive-current injection gain
        Iql1,   ///< Minimum reactive-current injection limit
        Iqh1,   ///< Maximum reactive-current injection limit
        Iqfrz,  ///< Post-dip frozen reactive-current injection
        Thld,   ///< Reactive-current injection hold time after dip recovery
        Thld2,  ///< Active-current limit hold time after dip recovery
        Qmax,   ///< Maximum reactive-power control limit
        Qmin,   ///< Minimum reactive-power control limit
        Kqp,    ///< Reactive-power proportional gain
        Kqi,    ///< Reactive-power integral gain
        Vmax,   ///< Maximum voltage-control limit
        Vmin,   ///< Minimum voltage-control limit
        Vref1,  ///< Voltage-reference bias
        Kvp,    ///< Voltage-control proportional gain
        Kvi,    ///< Voltage-control integral gain
        Tiq,    ///< Reactive-current command lag time constant
        Tpord,  ///< Active-power order filter time constant
        dPmax,  ///< Positive active-power ramp-rate limit
        dPmin,  ///< Negative active-power ramp-rate limit
        Pmax,   ///< Maximum active-power order limit
        Pmin,   ///< Minimum active-power order limit
        Imax,   ///< Maximum converter current
        Vq1,    ///< Reactive-current VDL voltage breakpoint 1
        Iq1,    ///< Reactive-current VDL limit 1
        Vq2,    ///< Reactive-current VDL voltage breakpoint 2
        Iq2,    ///< Reactive-current VDL limit 2
        Vq3,    ///< Reactive-current VDL voltage breakpoint 3
        Iq3,    ///< Reactive-current VDL limit 3
        Vq4,    ///< Reactive-current VDL voltage breakpoint 4
        Iq4,    ///< Reactive-current VDL limit 4
        Vp1,    ///< Active-current VDL voltage breakpoint 1
        Ip1,    ///< Active-current VDL limit 1
        Vp2,    ///< Active-current VDL voltage breakpoint 2
        Ip2,    ///< Active-current VDL limit 2
        Vp3,    ///< Active-current VDL voltage breakpoint 3
        Ip3,    ///< Active-current VDL limit 3
        Vp4,    ///< Active-current VDL voltage breakpoint 4
        Ip4     ///< Active-current VDL limit 4
      };

      /// Buses for the REECA electrical-control model.
      enum class ReecaBuses : size_t
      {
        bus, ///< Terminal bus ID
        SIZE
      };

      /// Signal inputs for the REECA electrical-control model.
      enum class ReecaSignalInputs : size_t
      {
        pe,     ///< Electrical active-power signal ID
        qgen,   ///< Reactive-power signal ID
        omegag, ///< Internal speed proxy signal ID; required when PFlag is 1
        qext,   ///< Optional reactive-power command signal ID
        pfaref, ///< Optional power-factor angle reference signal ID
        pref,   ///< Optional active-power reference signal ID
        SIZE
      };

      /// Signal outputs for the REECA electrical-control model.
      enum class ReecaSignalOutputs : size_t
      {
        iqcmd, ///< Reactive-current command output signal ID
        ipcmd, ///< Active-current command output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class ReecaMonitorableVariables
      {
        iqcmd, ///< Reactive-current command output
        ipcmd, ///< Active-current command output
        vmeas, ///< Filtered terminal voltage
        pmeas  ///< Filtered electrical power
      };

      template <typename real_type, typename index_type>
      struct ReecaData : public ComponentData<real_type,
                                              index_type,
                                              ReecaParameters,
                                              ReecaBuses,
                                              ReecaSignalInputs,
                                              ReecaSignalOutputs,
                                              ReecaMonitorableVariables>
      {
        ReecaData() = default;

        using Parameters           = ReecaParameters;
        using Buses                = ReecaBuses;
        using SignalInputs         = ReecaSignalInputs;
        using SignalOutputs        = ReecaSignalOutputs;
        using MonitorableVariables = ReecaMonitorableVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
