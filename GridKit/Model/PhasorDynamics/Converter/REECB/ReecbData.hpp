/**
 * @file ReecbData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REECB electrical-control model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Parameter keys for the REECB electrical-control model.
      enum class ReecbParameters
      {
        mva,    ///< Model MVA base
        PfFlag, ///< Power-factor control flag
        VFlag,  ///< Voltage-control mode flag
        QFlag,  ///< Reactive-power control flag
        Pqflag, ///< P/Q current-priority flag
        Trv,    ///< Voltage-measurement filter time constant
        Tp,     ///< Electrical-power measurement filter time constant
        Vref0,  ///< Outer-loop voltage reference
        Vdip,   ///< Low-voltage dip threshold
        Vup,    ///< High-voltage threshold
        dbd1,   ///< Overvoltage deadband threshold
        dbd2,   ///< Undervoltage deadband threshold
        kqv,    ///< Reactive-current injection gain
        Iql1,   ///< Minimum reactive-current injection limit
        Iqh1,   ///< Maximum reactive-current injection limit
        Qmax,   ///< Maximum reactive-power control limit
        Qmin,   ///< Minimum reactive-power control limit
        Kqp,    ///< Reactive-power proportional gain
        Kqi,    ///< Reactive-power integral gain
        Vmax,   ///< Maximum voltage-control limit
        Vmin,   ///< Minimum voltage-control limit
        Kvp,    ///< Voltage-control proportional gain
        Kvi,    ///< Voltage-control integral gain
        Tiq,    ///< Reactive-current command lag time constant
        Tpord,  ///< Active-power order filter time constant
        dPmax,  ///< Positive active-power ramp-rate limit
        dPmin,  ///< Negative active-power ramp-rate limit
        Pmax,   ///< Maximum active-power order limit
        Pmin,   ///< Minimum active-power order limit
        Imax    ///< Maximum converter current
      };

      /// Ports for the REECB electrical-control model.
      enum class ReecbPorts
      {
        bus,    ///< Terminal bus ID
        pe,     ///< Electrical active-power signal ID
        qgen,   ///< Reactive-power signal ID
        qext,   ///< Optional reactive-power command signal ID
        pfaref, ///< Optional power-factor angle reference signal ID
        pref,   ///< Optional active-power reference signal ID
        iqcmd,  ///< Reactive-current command output signal ID
        ipcmd   ///< Active-current command output signal ID
      };

      /// Variables available through the monitor interface.
      enum class ReecbMonitorableVariables
      {
        iqcmd,  ///< Reactive-current command output
        ipcmd,  ///< Active-current command output
        vmeas,  ///< Filtered terminal voltage
        pmeas,  ///< Filtered electrical power
        piq,    ///< Reactive-power PI controller state
        piv,    ///< Voltage PI controller state
        qv,     ///< Reactive-current command lag state
        pord,   ///< Filtered active-power order
        qref,   ///< Selected reactive-power reference
        sdip,   ///< Voltage-dip/overvoltage freeze indicator
        iqmax,  ///< Final reactive-current upper limit
        ipmax,  ///< Final active-current upper limit
        iqv,    ///< Reactive-current injection candidate
        vqctrl, ///< Reactive-power control PI output
        iqbase  ///< Base reactive-current command
      };

      template <typename RealT, typename IdxT>
      struct ReecbData : public ComponentData<RealT,
                                              IdxT,
                                              ReecbParameters,
                                              ReecbPorts,
                                              ReecbMonitorableVariables>
      {
        ReecbData() = default;

        using Parameters           = ReecbParameters;
        using Ports                = ReecbPorts;
        using MonitorableVariables = ReecbMonitorableVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
