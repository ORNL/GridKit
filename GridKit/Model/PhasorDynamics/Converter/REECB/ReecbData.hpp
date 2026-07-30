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
      /// Parameter keys for the REECB electrical-control model. `mva` is
      /// required; every other key is optional and retains its documented default.
      enum class ReecbParameters
      {
        mva,    ///< \f$S^\mathrm{base}\f$ REECB component power base
        PfFlag, ///< \f$s_\mathrm{pf}\f$ Power-factor control flag
        VFlag,  ///< \f$s_V\f$ Voltage-control mode flag
        QFlag,  ///< \f$s_Q\f$ Reactive-power control flag
        Pqflag, ///< \f$s_{PQ}\f$ P/Q current-priority flag
        Trv,    ///< \f$T_\mathrm{rv}\f$ Voltage-measurement time constant
        Tp,     ///< \f$T_\mathrm{p}\f$ Active-power measurement time constant
        Vref0,  ///< \f$V_0^\mathrm{ref}\f$ Outer-loop voltage reference
        Vdip,   ///< \f$V_\mathrm{dip}\f$ Low-voltage threshold
        Vup,    ///< \f$V_\mathrm{up}\f$ High-voltage threshold
        dbd1,   ///< \f$D_1^\mathrm{db}\f$ Lower voltage-error deadband
        dbd2,   ///< \f$D_2^\mathrm{db}\f$ Upper voltage-error deadband
        kqv,    ///< \f$K_\mathrm{qv}\f$ Reactive-current injection gain
        Iql1,   ///< \f$I_{q,\mathrm{inj}}^\min\f$ Minimum injection current
        Iqh1,   ///< \f$I_{q,\mathrm{inj}}^\max\f$ Maximum injection current
        Qmax,   ///< \f$Q^\max\f$ Maximum reactive-power control limit
        Qmin,   ///< \f$Q^\min\f$ Minimum reactive-power control limit
        Kqp,    ///< \f$K_\mathrm{qp}\f$ Reactive-power proportional gain
        Kqi,    ///< \f$K_\mathrm{qi}\f$ Reactive-power integral gain
        Vmax,   ///< \f$V^\max\f$ Maximum voltage-control limit
        Vmin,   ///< \f$V^\min\f$ Minimum voltage-control limit
        Kvp,    ///< \f$K_\mathrm{vp}\f$ Voltage-control proportional gain
        Kvi,    ///< \f$K_\mathrm{vi}\f$ Voltage-control integral gain
        Tiq,    ///< \f$T_\mathrm{iq}\f$ Reactive-current command time constant
        Tpord,  ///< \f$T_\mathrm{pord}\f$ Active-power order time constant
        dPmax,  ///< \f$R_P^\max\f$ Positive active-power ramp-rate limit
        dPmin,  ///< \f$R_P^\min\f$ Negative active-power ramp-rate limit
        Pmax,   ///< \f$P^\max\f$ Maximum active-power order limit
        Pmin,   ///< \f$P^\min\f$ Minimum active-power order limit
        Imax    ///< \f$I^\max\f$ Maximum converter current
      };

      /// Buses for the REECB electrical-control model.
      enum class ReecbBuses : size_t
      {
        bus, ///< Terminal bus ID
        SIZE
      };

      /// Optional signal inputs for the REECB electrical-control model.
      enum class ReecbSignalInputs : size_t
      {
        pe,     ///< Active-power feedback signal ID on system base
        qgen,   ///< Reactive-power feedback signal ID on system base
        qext,   ///< Reactive-power command signal ID on system base
        pfaref, ///< Power-factor angle reference signal ID in radians
        pref,   ///< Active-power reference signal ID on system base
        SIZE
      };

      /// Optional signal outputs for the REECB electrical-control model.
      enum class ReecbSignalOutputs : size_t
      {
        iqcmd, ///< Reactive-current command signal ID on system base
        ipcmd, ///< Active-current command signal ID on system base
        SIZE
      };

      /// Variables available through the monitor interface.
      enum class ReecbMonitorableVariables
      {
        iqcmd, ///< Reactive-current command on system base
        ipcmd, ///< Active-current command on system base
        vmeas, ///< Filtered terminal voltage
        pmeas  ///< Filtered active power on component base
      };

      /**
       * @brief Model data for the REECB controller: parameter values, the
       *        terminal bus, optional signals, and monitored-variable selections.
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
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
