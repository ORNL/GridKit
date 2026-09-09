/**
 * @file RepcaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for the REPCA plant-control model.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /// Parameter keys for `Repca`.
      enum class RepcaParameters
      {
        S,         ///< \f$S\f$ Rated apparent power [VA]
        V,         ///< \f$V\f$ Rated line-to-line RMS voltage [V]
        VcompFlag, ///< \f$s_\mathrm{comp}\f$ Voltage-compensation selector [boolean]
        RefFlag,   ///< \f$s_\mathrm{ref}\f$ Reactive-loop reference selector [boolean]
        Freqflag,  ///< \f$s_\mathrm{freq}\f$ Active-power output selector [boolean]
        Tfltr,     ///< \f$T_\mathrm{fltr}\f$ Voltage and reactive-power filter time constant [s]
        Vfrz,      ///< \f$V^\mathrm{frz}\f$ Reactive-power PI freeze threshold [p.u.]
        Rc,        ///< \f$R_c\f$ Line-drop resistance on component base [p.u.]
        Xc,        ///< \f$X_c\f$ Line-drop reactance on component base [p.u.]
        Kc,        ///< \f$K_c\f$ Reactive-droop coefficient [p.u.]
        dbdlow,    ///< \f$D_\mathrm{bd1}\f$ Lower reactive-loop deadband threshold [p.u.]
        dbdupper,  ///< \f$D_\mathrm{bd2}\f$ Upper reactive-loop deadband threshold [p.u.]
        emax,      ///< \f$e^{\max}\f$ Maximum reactive-loop error [p.u.]
        emin,      ///< \f$e^{\min}\f$ Minimum reactive-loop error [p.u.]
        Kp,        ///< \f$K_\mathrm{p}\f$ Reactive-power proportional gain [p.u.]
        Ki,        ///< \f$K_\mathrm{i}\f$ Reactive-power integral gain [p.u./s]
        Qmax,      ///< \f$Q^{\max}\f$ Maximum reactive-power command on component base [p.u.]
        Qmin,      ///< \f$Q^{\min}\f$ Minimum reactive-power command on component base [p.u.]
        Tft,       ///< \f$T_\mathrm{ft}\f$ Reactive-command lead time constant [s]
        Tfv,       ///< \f$T_\mathrm{fv}\f$ Reactive-command lag time constant [s]
        Tp,        ///< \f$T_\mathrm{p}\f$ Active-power measurement filter time constant [s]
        fdbd1,     ///< \f$D_\mathrm{bd1}^{f}\f$ Lower frequency-error deadband threshold [p.u.]
        fdbd2,     ///< \f$D_\mathrm{bd2}^{f}\f$ Upper frequency-error deadband threshold [p.u.]
        Ddn,       ///< \f$D_\mathrm{dn}\f$ Down-regulation (overfrequency) gain [p.u./p.u.]
        Dup,       ///< \f$D_\mathrm{up}\f$ Up-regulation (underfrequency) gain [p.u./p.u.]
        femax,     ///< \f$e_P^{\max}\f$ Maximum active-power error [p.u.]
        femin,     ///< \f$e_P^{\min}\f$ Minimum active-power error [p.u.]
        Kpg,       ///< \f$K_\mathrm{pg}\f$ Active-power proportional gain [p.u.]
        Kig,       ///< \f$K_\mathrm{ig}\f$ Active-power integral gain [p.u./s]
        Pmax,      ///< \f$P^{\max}\f$ Maximum active-power command on component base [p.u.]
        Pmin,      ///< \f$P^{\min}\f$ Minimum active-power command on component base [p.u.]
        Tlag       ///< \f$T_\mathrm{lag}\f$ Active-power command lag time constant [s]
      };

      enum class RepcaInputs : size_t
      {
        vd,      ///< \f$v_d\f$ Terminal d-axis voltage [V]
        vq,      ///< \f$v_q\f$ Terminal q-axis voltage [V]
        id,      ///< \f$i_d\f$ Terminal d-axis current [A]
        iq,      ///< \f$i_q\f$ Terminal q-axis current [A]
        freq,    ///< \f$f\f$ Absolute frequency [p.u.]
        vref,    ///< \f$V^{\mathrm{ref}}\f$ Voltage reference [V]
        pref,    ///< \f$P_{\mathrm{plant}}^{\mathrm{ref}}\f$ Active-power reference [W]
        qref,    ///< \f$Q^{\mathrm{ref}}\f$ Reactive-power reference [var]
        freqref, ///< \f$f^{\mathrm{ref}}\f$ Absolute frequency reference [p.u.]
        SIZE,
      };

      /// Signal outputs for `Repca`.
      enum class RepcaOutputs : size_t
      {
        qext, ///< \f$Q^\mathrm{ext}\f$ Reactive-power command [var]
        pext, ///< \f$P^\mathrm{ext}\f$ Active-power command [W]
        SIZE  ///< Number of REPCA signal-output ports
      };

      /// Variables available through the monitor interface.
      enum class RepcaMonitorableVariables
      {
        qext,  ///< \f$Q^\mathrm{ext}\f$ Reactive-power command output [var]
        pext,  ///< \f$P^\mathrm{ext}\f$ Active-power command output [W]
        vmeas, ///< \f$V^\mathrm{meas}\f$ Filtered regulated voltage [p.u.]
        qmeas, ///< \f$Q^\mathrm{meas}\f$ Filtered reactive-power signal on component base [p.u.]
        pmeas  ///< \f$P^\mathrm{meas}\f$ Filtered active-power signal on component base [p.u.]
      };

      /**
       * @brief Model data for REPCA parameters, signal ports, and monitors.
       *
       * @tparam real_type Real parameter value type.
       * @tparam index_type Integer index type.
       *
       * @see Repca
       */
      template <typename real_type, typename index_type>
      struct RepcaData : public ComponentData<real_type,
                                              index_type,
                                              RepcaParameters,
                                              RepcaInputs,
                                              RepcaOutputs,
                                              RepcaMonitorableVariables>
      {
        RepcaData() = default;

        using Parameters           = RepcaParameters;
        using Inputs               = RepcaInputs;
        using Outputs              = RepcaOutputs;
        using MonitorableVariables = RepcaMonitorableVariables;
      };
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
