/**
 * @file ReecbData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REECB electrical-control model.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Parameters for REECB.
      enum class ReecbParameters
      {
        S,      ///< \f$S\f$ Rated apparent power [VA]
        V,      ///< \f$V\f$ Rated line-to-line RMS voltage [V]
        PfFlag, ///< \f$s_\mathrm{pf}\f$ Power-factor control selector: true = power-factor control, false = reactive-power control [boolean]
        VFlag,  ///< \f$s_V\f$ Voltage-reference selector under \f$s_Q=1\f$: true = cascaded Q-PI voltage command, false = direct external voltage reference [boolean]
        QFlag,  ///< \f$s_Q\f$ Reactive-path selector: true = Volt/VAr PI control, false = reactive-current lag [boolean]
        Pqflag, ///< \f$s_\mathrm{pq}\f$ Converter current-priority selector: true = P priority, false = Q priority [boolean]
        Trv,    ///< \f$T_\mathrm{rv}\f$ Voltage-measurement filter time constant [s]
        Tp,     ///< \f$T_\mathrm{p}\f$ Electrical-power measurement filter time constant [s]
        Vref0,  ///< \f$V_0^{\mathrm{ref}}\f$ Reactive-current-injection voltage reference [p.u.]
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
        Tiq,    ///< \f$T_\mathrm{iq}\f$ Reactive-current command lag time constant [s]
        Tpord,  ///< \f$T_\mathrm{pord}\f$ Active-power order filter time constant [s]
        dPmax,  ///< \f$R_P^\max\f$ Positive active-power ramp-rate limit on component base [p.u./s]
        dPmin,  ///< \f$R_P^\min\f$ Negative active-power ramp-rate limit on component base [p.u./s]
        Pmax,   ///< \f$P^\max\f$ Maximum active-power order limit on component base [p.u.]
        Pmin,   ///< \f$P^\min\f$ Minimum active-power order limit on component base [p.u.]
        Imax    ///< \f$I^\max\f$ Maximum converter current on component base [p.u.]
      };

      enum class ReecbInputs : size_t
      {
        vd,     ///< \f$v_d\f$ Terminal d-axis voltage [V]
        vq,     ///< \f$v_q\f$ Terminal q-axis voltage [V]
        id,     ///< \f$i_d\f$ Terminal d-axis current [A]
        iq,     ///< \f$i_q\f$ Terminal q-axis current [A]
        Pref,   ///< \f$P^{\mathrm{ref}}\f$ Active-power reference [W]
        Qref,   ///< \f$Q^{\mathrm{ref}}\f$ Reactive-power reference [var]
        Vref,   ///< \f$V^{\mathrm{ref}}\f$ Terminal-voltage reference [V]
        pfaref, ///< \f$\phi^{\mathrm{ref}}\f$ Power-factor angle reference [rad]
        SIZE
      };

      enum class ReecbOutputs : size_t
      {
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Terminal d-axis current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Terminal q-axis current command [A]
        SIZE
      };

      enum class ReecbMonitorableVariables
      {
        icmdd, ///< \f$i_d^{\mathrm{cmd}}\f$ Terminal d-axis current command [A]
        icmdq, ///< \f$i_q^{\mathrm{cmd}}\f$ Terminal q-axis current command [A]
        ipcmd, ///< \f$I_p^{\mathrm{cmd}}\f$ Active-current command on component base [p.u.]
        iqcmd, ///< \f$I_q^{\mathrm{cmd}}\f$ Reactive-current command on component base [p.u.]
        iqv,   ///< \f$I_q^{\mathrm{inj}}\f$ Supplementary reactive current on component base [p.u.]
        vmeas, ///< \f$V^{\mathrm{meas}}\f$ Filtered terminal voltage [p.u.]
        pmeas  ///< \f$P^{\mathrm{meas}}\f$ Filtered active power on component base [p.u.]
      };

      /**
       * @brief Model data for REECB parameters, signal ports, and monitored variables.
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
                                              ReecbInputs,
                                              ReecbOutputs,
                                              ReecbMonitorableVariables>
      {
        ReecbData() = default;

        using Parameters           = ReecbParameters;
        using Inputs               = ReecbInputs;
        using Outputs              = ReecbOutputs;
        using MonitorableVariables = ReecbMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
