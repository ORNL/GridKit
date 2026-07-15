/**
 * @file RepcaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REPCA converter plant-control model.
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
        mva,       ///< REPCA component power base
        VcompFlag, ///< Voltage-compensation mode flag
        RefFlag,   ///< Reactive-loop reference flag
        Freqflag,  ///< Active-power control output flag
        Tfltr,     ///< Voltage and reactive-power measurement filter time constant
        Vfrz,      ///< Reactive-PI freeze voltage threshold
        Rc,        ///< Line-drop compensation resistance
        Xc,        ///< Line-drop compensation reactance
        Kc,        ///< Reactive-current compensation gain
        dbdlow,    ///< Lower reactive-loop deadband threshold
        dbdupper,  ///< Upper reactive-loop deadband threshold
        emax,      ///< Maximum reactive-loop error limit
        emin,      ///< Minimum reactive-loop error limit
        Kp,        ///< Reactive-power controller proportional gain
        Ki,        ///< Reactive-power controller integral gain
        Qmax,      ///< Maximum reactive-power command
        Qmin,      ///< Minimum reactive-power command
        Tft,       ///< Reactive-command lead time constant
        Tfv,       ///< Reactive-command lag time constant
        Tp,        ///< Active-power measurement filter time constant
        fdbd1,     ///< Lower frequency-error deadband threshold
        fdbd2,     ///< Upper frequency-error deadband threshold
        Ddn,       ///< Down-frequency droop response gain
        Dup,       ///< Up-frequency droop response gain
        femax,     ///< Maximum active-power error limit
        femin,     ///< Minimum active-power error limit
        Kpg,       ///< Active-power controller proportional gain
        Kig,       ///< Active-power controller integral gain
        Pmax,      ///< Maximum active-power command
        Pmin,      ///< Minimum active-power command
        Tlag       ///< Active-power command lag time constant
      };

      /// Buses for the REPCA plant-control model.
      enum class RepcaBuses : size_t
      {
        bus, ///< Regulated bus ID
        SIZE
      };

      /// Signal inputs for the REPCA plant-control model.
      enum class RepcaSignalInputs : size_t
      {
        ibranchr,  ///< Branch current real-component signal ID
        ibranchi,  ///< Branch current imaginary-component signal ID
        pbranch,   ///< Branch active-power signal ID
        qbranch,   ///< Branch reactive-power signal ID
        freq,      ///< Frequency input signal ID
        freqref,   ///< Frequency-reference signal ID
        vref,      ///< Voltage-reference signal ID
        qref,      ///< Reactive-power reference signal ID
        pplantref, ///< Plant active-power reference signal ID
        SIZE
      };

      /// Signal outputs for the REPCA plant-control model.
      enum class RepcaSignalOutputs : size_t
      {
        qext, ///< Reactive-power command output signal ID
        pext, ///< Active-power command output signal ID
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class RepcaMonitorableVariables
      {
        qext,  ///< Reactive-power command output
        pext,  ///< Active-power command output
        vmeas, ///< Filtered regulated voltage
        qmeas, ///< Filtered reactive-power signal
        pmeas  ///< Filtered active-power signal
      };

      template <typename real_type, typename index_type>
      struct RepcaData : public ComponentData<real_type,
                                              index_type,
                                              RepcaParameters,
                                              RepcaBuses,
                                              RepcaSignalInputs,
                                              RepcaSignalOutputs,
                                              RepcaMonitorableVariables>
      {
        RepcaData() = default;

        using Parameters           = RepcaParameters;
        using Buses                = RepcaBuses;
        using SignalInputs         = RepcaSignalInputs;
        using SignalOutputs        = RepcaSignalOutputs;
        using MonitorableVariables = RepcaMonitorableVariables;
      };
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
