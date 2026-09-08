#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class RegfmaParameters
    {
      S,      ///< \f$S_\mathrm{b}\f$ Rated three-phase apparent power [VA]
      V,      ///< \f$V_\mathrm{b}\f$ Line-to-line RMS voltage base [V]
      omega0, ///< \f$\omega_0\f$ Rated angular frequency [rad/s]
      XL,     ///< \f$X_L\f$ Coupling reactance [p.u.]
      mp,     ///< \f$m_p\f$ Active-power droop [p.u.]
      mq,     ///< \f$m_q\f$ Reactive-power droop [p.u.]
      kpv,    ///< \f$k_{\mathrm{pv}}\f$ Voltage proportional gain [p.u.]
      kiv,    ///< \f$k_{\mathrm{iv}}\f$ Voltage integral gain [\f$\mathrm{s}^{-1}\f$]
      Emin,   ///< \f$E^{\min}\f$ Minimum internal voltage [p.u.]
      Emax,   ///< \f$E^{\max}\f$ Maximum internal voltage [p.u.]
      Pmin,   ///< \f$P^{\min}\f$ Minimum active power [p.u.]
      Pmax,   ///< \f$P^{\max}\f$ Maximum active power [p.u.]
      Qmin,   ///< \f$Q^{\min}\f$ Minimum reactive power [p.u.]
      Qmax,   ///< \f$Q^{\max}\f$ Maximum reactive power [p.u.]
      kppmax, ///< \f$k_{\mathrm{ppmax}}\f$ Active-power limit proportional gain [p.u.]
      kipmax, ///< \f$k_{\mathrm{ipmax}}\f$ Active-power limit integral gain [\f$\mathrm{s}^{-1}\f$]
      kpqmax, ///< \f$k_{\mathrm{pqmax}}\f$ Reactive-power limit proportional gain [p.u.]
      kiqmax, ///< \f$k_{\mathrm{iqmax}}\f$ Reactive-power limit integral gain [\f$\mathrm{s}^{-1}\f$]
      TPf,    ///< \f$T_{Pf}\f$ Active-power measurement time constant [s]
      TQf,    ///< \f$T_{Qf}\f$ Reactive-power measurement time constant [s]
      TVf,    ///< \f$T_{Vf}\f$ Voltage measurement time constant [s]
      ImaxF,  ///< \f$I_F^{\max}\f$ Maximum transient current [p.u.]
      VFlag,  ///< \f$\mathrm{VFlag}\f$ Regulate terminal voltage when true, internal voltage otherwise
      QVFlag, ///< \f$\mathrm{QVFlag}\f$ Infer the voltage reference when true, reactive reference otherwise
    };

    enum class RegfmaInputs : size_t
    {
      va,   ///< \f$v_a\f$ Phase-a terminal voltage [V]
      vb,   ///< \f$v_b\f$ Phase-b terminal voltage [V]
      vc,   ///< \f$v_c\f$ Phase-c terminal voltage [V]
      pref, ///< \f$P_\mathrm{ref}\f$ Active-power reference [p.u.]
      qref, ///< \f$Q_\mathrm{ref}\f$ Reactive-power reference [p.u.]
      vref, ///< \f$V_\mathrm{ref}\f$ Voltage reference [p.u.]
      SIZE,
    };

    enum class RegfmaOutputs : size_t
    {
      ia, ///< \f$i_a\f$ Phase-a current injection [A]
      ib, ///< \f$i_b\f$ Phase-b current injection [A]
      ic, ///< \f$i_c\f$ Phase-c current injection [A]
      SIZE,
    };

    enum class RegfmaMonitorableVariables
    {
      pf,     ///< \f$P_f\f$ Filtered active power [p.u.]
      qf,     ///< \f$Q_f\f$ Filtered reactive power [p.u.]
      vf,     ///< \f$V_f\f$ Filtered voltage magnitude [p.u.]
      xpmax,  ///< \f$x_P^{\max}\f$ Upper active-power-limit integral [p.u.]
      xpmin,  ///< \f$x_P^{\min}\f$ Lower active-power-limit integral [p.u.]
      xqmax,  ///< \f$x_Q^{\max}\f$ Upper reactive-power-limit integral [p.u.]
      xqmin,  ///< \f$x_Q^{\min}\f$ Lower reactive-power-limit integral [p.u.]
      xv,     ///< \f$x_V\f$ Voltage-control integral [p.u.]
      delta,  ///< \f$\delta\f$ Internal angle deviation [rad]
      ia,     ///< \f$i_a\f$ Phase-a current injection [A]
      ib,     ///< \f$i_b\f$ Phase-b current injection [A]
      ic,     ///< \f$i_c\f$ Phase-c current injection [A]
      ea,     ///< \f$e_a\f$ Corrected phase-a source voltage [V]
      eb,     ///< \f$e_b\f$ Corrected phase-b source voltage [V]
      ec,     ///< \f$e_c\f$ Corrected phase-c source voltage [V]
      omega,  ///< \f$\omega\f$ Internal angular frequency [rad/s]
      edroop, ///< \f$E_\mathrm{droop}\f$ Droop voltage magnitude [p.u.]
      p,      ///< \f$S_\mathrm{b}P\f$ Terminal active power [W]
      q,      ///< \f$S_\mathrm{b}Q\f$ Terminal reactive power [var]
      v,      ///< \f$V_\mathrm{b}V\f$ Terminal voltage magnitude [V]
    };

    template <typename real_type, typename index_type>
    struct RegfmaData : public ComponentData<real_type,
                                             index_type,
                                             RegfmaParameters,
                                             RegfmaInputs,
                                             RegfmaOutputs,
                                             RegfmaMonitorableVariables>
    {
      RegfmaData() = default;

      using Parameters           = RegfmaParameters;
      using Inputs               = RegfmaInputs;
      using Outputs              = RegfmaOutputs;
      using MonitorableVariables = RegfmaMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
