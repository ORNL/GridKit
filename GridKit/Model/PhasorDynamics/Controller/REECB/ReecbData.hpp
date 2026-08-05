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
    namespace Controller
    {
      /// Parameters for REECB.
      enum class ReecbParameters
      {
        mva,    ///< \f$S^\mathrm{base}\f$ Component power base [MVA]
        PfFlag, ///< \f$s_\mathrm{pf}\f$ Power-factor control selector: true = power-factor control, false = reactive-power control [boolean]
        VFlag,  ///< \f$s_V\f$ Voltage-reference selector under \f$s_Q=1\f$: true = cascaded Q-PI voltage command, false = direct external voltage reference [boolean]
        QFlag,  ///< \f$s_Q\f$ Reactive-path selector: true = Volt/VAr PI control, false = reactive-current lag [boolean]
        Pqflag, ///< \f$s_{PQ}\f$ Converter current-priority selector: true = P priority, false = Q priority [boolean]
        Trv,    ///< \f$T_\mathrm{rv}\f$ Voltage-measurement filter time constant [sec]
        Tp,     ///< \f$T_\mathrm{p}\f$ Electrical-power measurement filter time constant [sec]
        Vref0,  ///< \f$V^\mathrm{ref}\f$ Reactive-current-injection voltage reference [p.u.]
        Vdip,   ///< \f$V_\mathrm{dip}\f$ Low-voltage threshold for the voltage-band gate [p.u.]
        Vup,    ///< \f$V_\mathrm{up}\f$ High-voltage threshold for the voltage-band gate [p.u.]
        dbd1,   ///< \f$D_1^\mathrm{db}\f$ Lower voltage-error deadband threshold [p.u.]
        dbd2,   ///< \f$D_2^\mathrm{db}\f$ Upper voltage-error deadband threshold [p.u.]
        kqv,    ///< \f$K_\mathrm{qv}\f$ Reactive-current injection gain [p.u.]
        Iql1,   ///< \f$I_{q,\mathrm{inj}}^\min\f$ Minimum reactive-current injection on component base [p.u.]
        Iqh1,   ///< \f$I_{q,\mathrm{inj}}^\max\f$ Maximum reactive-current injection on component base [p.u.]
        Qmax,   ///< \f$Q^\max\f$ Maximum reactive-power control output on component base [p.u.]
        Qmin,   ///< \f$Q^\min\f$ Minimum reactive-power control output on component base [p.u.]
        Kqp,    ///< \f$K_\mathrm{qp}\f$ Reactive-power proportional gain [p.u.]
        Kqi,    ///< \f$K_\mathrm{qi}\f$ Reactive-power integral gain [p.u./s]
        Vmax,   ///< \f$V^\max\f$ Maximum voltage-control output [p.u.]
        Vmin,   ///< \f$V^\min\f$ Minimum voltage-control output [p.u.]
        Kvp,    ///< \f$K_\mathrm{vp}\f$ Voltage-control proportional gain [p.u.]
        Kvi,    ///< \f$K_\mathrm{vi}\f$ Voltage-control integral gain [p.u./s]
        Tiq,    ///< \f$T_\mathrm{iq}\f$ Reactive-current command lag time constant [sec]
        Tpord,  ///< \f$T_\mathrm{pord}\f$ Active-power order filter time constant [sec]
        dPmax,  ///< \f$R_P^\max\f$ Positive active-power ramp-rate limit on component base [p.u./s]
        dPmin,  ///< \f$R_P^\min\f$ Negative active-power ramp-rate limit on component base [p.u./s]
        Pmax,   ///< \f$P^\max\f$ Maximum active-power order limit on component base [p.u.]
        Pmin,   ///< \f$P^\min\f$ Minimum active-power order limit on component base [p.u.]
        Imax    ///< \f$I^\max\f$ Maximum converter current on component base [p.u.]
      };

      /// Buses for the REECB electrical-control model.
      enum class ReecbBuses : size_t
      {
        bus, ///< \f$V_\mathrm{r},V_\mathrm{i}\f$ Required Known terminal-bus voltage [p.u.]
        SIZE ///< Number of REECB bus ports
      };

      /// Signal inputs for the REECB electrical-control model.
      enum class ReecbSignalInputs : size_t
      {
        pe,     ///< \f$P_e\f$ Optional Known active-power feedback input on system base [p.u.]
        qgen,   ///< \f$Q^\mathrm{gen}\f$ Optional Known reactive-power feedback input on system base [p.u.]
        qext,   ///< \f$Q^\mathrm{ext}\f$ Optional Unknown Volt/VAr reference input: system-base reactive power [p.u.], or the terminal-voltage reference [p.u.] when \f$s_Q=1\f$ and \f$s_V=0\f$
        pfaref, ///< \f$\phi^\mathrm{ref}\f$ Optional Unknown power-factor angle-reference input [rad]
        pref,   ///< \f$P^\mathrm{ref}\f$ Optional Unknown active-power reference input on system base [p.u.]
        SIZE    ///< Number of REECB signal-input ports
      };

      /// Signal outputs for the REECB electrical-control model.
      enum class ReecbSignalOutputs : size_t
      {
        iqcmd, ///< \f$I_q^\mathrm{cmd}\f$ Optional Known reactive-current command output on system base [p.u.]
        ipcmd, ///< \f$I_p^\mathrm{cmd}\f$ Optional Known active-current command output on system base [p.u.]
        SIZE   ///< Number of REECB signal-output ports
      };

      /// Variables available through the monitor interface.
      enum class ReecbMonitorableVariables
      {
        iqcmd, ///< \f$I_q^\mathrm{cmd}\f$ Reactive-current command output on system base [p.u.]
        ipcmd, ///< \f$I_p^\mathrm{cmd}\f$ Active-current command output on system base [p.u.]
        vmeas, ///< \f$V^\mathrm{meas}\f$ Filtered terminal voltage [p.u.]
        pmeas  ///< \f$P^\mathrm{meas}\f$ Filtered electrical power on component base [p.u.]
      };

      /**
       * @brief Model data for REECB parameters, bus and signal ports, and monitored variables.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see Reecb
       */
      template <typename real_type, typename index_type>
      struct ReecbData : public ComponentData<real_type,
                                              index_type,
                                              ReecbParameters,
                                              ReecbBuses,
                                              ReecbSignalInputs,
                                              ReecbSignalOutputs,
                                              ReecbMonitorableVariables>
      {
        ReecbData() = default;

        using Parameters           = ReecbParameters;
        using Buses                = ReecbBuses;
        using SignalInputs         = ReecbSignalInputs;
        using SignalOutputs        = ReecbSignalOutputs;
        using MonitorableVariables = ReecbMonitorableVariables;
      };
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
