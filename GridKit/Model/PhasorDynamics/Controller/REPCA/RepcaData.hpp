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
    namespace Controller
    {
      /// Parameter keys for `Repca`.
      enum class RepcaParameters
      {
        mva,       ///< \f$S^\mathrm{base}\f$ Component power base [MVA]
        VcompFlag, ///< \f$s_\mathrm{comp}\f$ Voltage-compensation selector [boolean]
        RefFlag,   ///< \f$s_\mathrm{ref}\f$ Reactive-loop reference selector [boolean]
        Freqflag,  ///< \f$s_\mathrm{freq}\f$ Active-power output selector [boolean]
        Tfltr,     ///< \f$T_\mathrm{fltr}\f$ Voltage and reactive-power filter time constant [sec]
        Vfrz,      ///< \f$V^\mathrm{frz}\f$ Reactive-power PI freeze threshold [p.u.]
        Rc,        ///< \f$R_c\f$ Line-drop resistance on component base [p.u.]
        Xc,        ///< \f$X_c\f$ Line-drop reactance on component base [p.u.]
        Kc,        ///< \f$K_c\f$ Reactive-current compensation gain [p.u.]
        dbdlow,    ///< \f$D_\mathrm{bd1}\f$ Lower reactive-loop deadband threshold [p.u.]
        dbdupper,  ///< \f$D_\mathrm{bd2}\f$ Upper reactive-loop deadband threshold [p.u.]
        emax,      ///< \f$e^{\max}\f$ Maximum reactive-loop error [p.u.]
        emin,      ///< \f$e^{\min}\f$ Minimum reactive-loop error [p.u.]
        Kp,        ///< \f$K_\mathrm{p}\f$ Reactive-power proportional gain [p.u.]
        Ki,        ///< \f$K_\mathrm{i}\f$ Reactive-power integral gain [p.u./s]
        Qmax,      ///< \f$Q^{\max}\f$ Maximum reactive-power command on component base [p.u.]
        Qmin,      ///< \f$Q^{\min}\f$ Minimum reactive-power command on component base [p.u.]
        Tft,       ///< \f$T_\mathrm{ft}\f$ Reactive-command lead time constant [sec]
        Tfv,       ///< \f$T_\mathrm{fv}\f$ Reactive-command lag time constant [sec]
        Tp,        ///< \f$T_\mathrm{p}\f$ Active-power measurement filter time constant [sec]
        fdbd1,     ///< \f$D_\mathrm{bd1}^{f}\f$ Lower frequency-error deadband threshold [p.u.]
        fdbd2,     ///< \f$D_\mathrm{bd2}^{f}\f$ Upper frequency-error deadband threshold [p.u.]
        Ddn,       ///< \f$D_\mathrm{dn}\f$ Down-frequency droop response gain [p.u.]
        Dup,       ///< \f$D_\mathrm{up}\f$ Up-frequency droop response gain [p.u.]
        femax,     ///< \f$e_P^{\max}\f$ Maximum active-power error [p.u.]
        femin,     ///< \f$e_P^{\min}\f$ Minimum active-power error [p.u.]
        Kpg,       ///< \f$K_\mathrm{pg}\f$ Active-power proportional gain [p.u.]
        Kig,       ///< \f$K_\mathrm{ig}\f$ Active-power integral gain [p.u./s]
        Pmax,      ///< \f$P^{\max}\f$ Maximum active-power command on component base [p.u.]
        Pmin,      ///< \f$P^{\min}\f$ Minimum active-power command on component base [p.u.]
        Tlag       ///< \f$T_\mathrm{lag}\f$ Active-power command lag time constant [sec]
      };

      /// Buses for `Repca`.
      enum class RepcaBuses : size_t
      {
        bus, ///< \f$V_\mathrm{r},V_\mathrm{i}\f$ Required Known regulated-bus voltage [p.u.]
        SIZE ///< Number of REPCA bus ports
      };

      /// Signal inputs for the `Repca`.
      enum class RepcaSignalInputs : size_t
      {
        ir,      ///< \f$I_\mathrm{r}\f$ Required Known branch-current real input on system base [p.u.]
        ii,      ///< \f$I_\mathrm{i}\f$ Required Known branch-current imaginary input on system base [p.u.]
        p,       ///< \f$P\f$ Required Known branch active-power input on system base [p.u.]
        q,       ///< \f$Q\f$ Required Known branch reactive-power input on system base [p.u.]
        freq,    ///< \f$f\f$ Optional Known frequency input [p.u.]
        vref,    ///< \f$V^\mathrm{ref}\f$ Optional Unknown voltage-reference input [p.u.]
        pref,    ///< \f$P_\mathrm{plant}^\mathrm{ref}\f$ Optional Unknown plant active-power reference on system base [p.u.]
        qref,    ///< \f$Q^\mathrm{ref}\f$ Optional Unknown reactive-power reference on system base [p.u.]
        freqref, ///< \f$f^\mathrm{ref}\f$ Optional Unknown frequency-reference input [p.u.]
        SIZE     ///< Number of REPCA signal-input ports
      };

      /// Signal outputs for `Repca`.
      enum class RepcaSignalOutputs : size_t
      {
        qext, ///< \f$Q^\mathrm{ext}\f$ Optional Known reactive-power command output on system base [p.u.]
        pext, ///< \f$P^\mathrm{ext}\f$ Optional Known active-power command output on system base [p.u.]
        SIZE  ///< Number of REPCA signal-output ports
      };

      /// Variables available through the monitor interface.
      enum class RepcaMonitorableVariables
      {
        qext,  ///< \f$Q^\mathrm{ext}\f$ Reactive-power command output on system base [p.u.]
        pext,  ///< \f$P^\mathrm{ext}\f$ Active-power command output on system base [p.u.]
        vmeas, ///< \f$V^\mathrm{meas}\f$ Filtered regulated voltage [p.u.]
        qmeas, ///< \f$Q^\mathrm{meas}\f$ Filtered reactive-power signal on component base [p.u.]
        pmeas  ///< \f$P^\mathrm{meas}\f$ Filtered active-power signal on component base [p.u.]
      };

      /**
       * @brief Model data for REPCA parameters, bus and signal ports, and monitored variables.
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
    } // namespace Controller
  } // namespace PhasorDynamics
} // namespace GridKit
