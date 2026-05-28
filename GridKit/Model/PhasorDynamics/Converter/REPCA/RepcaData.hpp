/**
 * @file RepcaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REPCA plant-control model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /// Parameter keys for the REPCA plant-control model.
      enum class RepcaParameters
      {
        mva,       ///< Component MVA base; zero uses system base
        VcompFlag, ///< Voltage-compensation mode flag
        RefFlag,   ///< Reactive-loop reference flag
        Freqflag,  ///< Active-power output flag
        Tfltr,     ///< Voltage and reactive-power filter time constant
        Tft,       ///< Reactive-command lead time constant
        Tfv,       ///< Reactive-command lag time constant
        Tp,        ///< Active-power measurement filter time constant
        Tlag,      ///< Active-power command lag time constant
        Vfrz,      ///< Reactive PI freeze voltage threshold
        Rc,        ///< Line-drop compensation resistance
        Xc,        ///< Line-drop compensation reactance
        Kc,        ///< Reactive-current compensation gain
        dbdlow,    ///< Lower reactive-loop deadband threshold
        dbdupper,  ///< Upper reactive-loop deadband threshold
        emax,      ///< Maximum reactive-loop error limit
        emin,      ///< Minimum reactive-loop error limit
        Kp,        ///< Reactive controller proportional gain
        Ki,        ///< Reactive controller integral gain
        Qmax,      ///< Maximum reactive-power command
        Qmin,      ///< Minimum reactive-power command
        fdbd1,     ///< Lower frequency-error deadband threshold
        fdbd2,     ///< Upper frequency-error deadband threshold
        Ddn,       ///< Down-frequency droop response gain
        Dup,       ///< Up-frequency droop response gain
        femax,     ///< Maximum active-power error limit
        femin,     ///< Minimum active-power error limit
        Kpg,       ///< Active-power proportional gain
        Kig,       ///< Active-power integral gain
        Pmax,      ///< Maximum active-power command
        Pmin       ///< Minimum active-power command
      };

      /// Ports for the REPCA plant-control model.
      enum class RepcaPorts
      {
        bus,       ///< Regulated bus ID
        ibranchr,  ///< Branch current real signal ID
        ibranchi,  ///< Branch current imaginary signal ID
        qbranch,   ///< Branch reactive-power signal ID
        pbranch,   ///< Branch active-power signal ID
        vref,      ///< Optional voltage reference signal ID
        qref,      ///< Optional reactive-power reference signal ID
        pplantref, ///< Optional plant active-power reference signal ID
        freq,      ///< Optional frequency input signal ID
        freqref,   ///< Optional frequency reference signal ID
        qext,      ///< Reactive-power command output signal ID
        pext       ///< Active-power command output signal ID
      };

      /// Variables available through the monitor interface.
      enum class RepcaMonitorableVariables
      {
        qext,  ///< Reactive-power command output
        pext,  ///< Active-power command output
        vmeas, ///< Filtered regulated voltage
        qmeas, ///< Filtered reactive-power signal
        pmeas, ///< Filtered active-power signal
        pref,  ///< Active-power command lag state
        vctrl, ///< Selected voltage-measurement input
        sfrz,  ///< Reactive PI voltage-enable indicator
        qpi,   ///< Reactive PI output
        pfreq, ///< Frequency droop active-power correction
        ppi    ///< Active-power PI output
      };

      template <typename real_type, typename index_type>
      struct RepcaData : public ComponentData<real_type,
                                              index_type,
                                              RepcaParameters,
                                              RepcaPorts,
                                              RepcaMonitorableVariables>
      {
        RepcaData() = default;

        using Parameters           = RepcaParameters;
        using Ports                = RepcaPorts;
        using MonitorableVariables = RepcaMonitorableVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
