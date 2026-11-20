/**
 * @file RepcaData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT renewable energy plant control.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Repca`
    enum class RepcaInputs : size_t
    {
      v,       ///< \f$V\f$ Regulated voltage, line-to-line RMS [V]
      p,       ///< \f$P\f$ Active power into the grid [W]
      q,       ///< \f$Q\f$ Reactive power into the grid [var]
      omega,   ///< \f$\omega\f$ Frequency [rad/s]
      vref,    ///< \f$V^\mathrm{ref}\f$ Voltage reference [p.u.]
      pref,    ///< \f$P^\mathrm{ref}\f$ Active-power reference on the system power base [p.u.]
      qref,    ///< \f$Q^\mathrm{ref}\f$ Reactive-power reference on the system power base [p.u.]
      freqref, ///< \f$f^\mathrm{ref}\f$ Frequency reference [p.u.]
    };

    /// Parameters of a `Repca`
    enum class RepcaParameters : size_t
    {
      mva,       ///< \f$S\f$ Three-phase device rating [MVA]
      VcompFlag, ///< \f$s_\mathrm{comp}\f$ Voltage-compensation selector
      RefFlag,   ///< \f$s_\mathrm{ref}\f$ Reactive-loop reference selector
      Freqflag,  ///< \f$s_\mathrm{freq}\f$ Active-power output selector
      Tfltr,     ///< \f$T_\mathrm{fltr}\f$ Voltage and reactive-power filter time constant [s]
      Vfrz,      ///< \f$V^\mathrm{frz}\f$ Reactive-power PI freeze-voltage threshold [p.u.]
      Rc,        ///< \f$R_c\f$ Line-drop compensation resistance [p.u.]
      Xc,        ///< \f$X_c\f$ Line-drop compensation reactance [p.u.]
      Kc,        ///< \f$K_c\f$ Reactive-current compensation gain [p.u.]
      dbdlow,    ///< \f$D_\mathrm{bd1}\f$ Lower reactive-loop deadband threshold [p.u.]
      dbdupper,  ///< \f$D_\mathrm{bd2}\f$ Upper reactive-loop deadband threshold [p.u.]
      emax,      ///< \f$e^{\max}\f$ Upper limit on the signed reactive-loop error [p.u.]
      emin,      ///< \f$e^{\min}\f$ Lower limit on the signed reactive-loop error [p.u.]
      Kp,        ///< \f$K_\mathrm{p}\f$ Reactive-power controller proportional gain [p.u.]
      Ki,        ///< \f$K_\mathrm{i}\f$ Reactive-power controller integral gain [p.u./s]
      Qmax,      ///< \f$Q^{\max}\f$ Maximum reactive-power command [p.u.]
      Qmin,      ///< \f$Q^{\min}\f$ Minimum reactive-power command [p.u.]
      Tft,       ///< \f$T_\mathrm{ft}\f$ Reactive-command lead time constant [s]
      Tfv,       ///< \f$T_\mathrm{fv}\f$ Reactive-command lag time constant [s]
      Tp,        ///< \f$T_\mathrm{p}\f$ Active-power measurement filter time constant [s]
      fdbd1,     ///< \f$D_\mathrm{bd1}^{f}\f$ Lower frequency-error deadband threshold [p.u.]
      fdbd2,     ///< \f$D_\mathrm{bd2}^{f}\f$ Upper frequency-error deadband threshold [p.u.]
      Ddn,       ///< \f$D_\mathrm{dn}\f$ Down-regulation gain for overfrequency [p.u.]
      Dup,       ///< \f$D_\mathrm{up}\f$ Up-regulation gain for underfrequency [p.u.]
      femax,     ///< \f$e_P^{\max}\f$ Upper limit on the signed active-power error [p.u.]
      femin,     ///< \f$e_P^{\min}\f$ Lower limit on the signed active-power error [p.u.]
      Kpg,       ///< \f$K_\mathrm{pg}\f$ Active-power controller proportional gain [p.u.]
      Kig,       ///< \f$K_\mathrm{ig}\f$ Active-power controller integral gain [p.u./s]
      Pmax,      ///< \f$P^{\max}\f$ Maximum active-power command [p.u.]
      Pmin,      ///< \f$P^{\min}\f$ Minimum active-power command [p.u.]
      Tlag,      ///< \f$T_\mathrm{lag}\f$ Active-power command lag time constant [s]
      V,         ///< \f$V\f$ Rated line-to-line RMS voltage [V]
    };

    /// Data of a `Repca`
    template <typename RealT, typename IdxT>
    struct RepcaData
    {
      using Parameters = RepcaParameters;

      RealT mva{100.0};
      bool  VcompFlag{true};
      bool  RefFlag{true};
      bool  Freqflag{false};
      RealT Tfltr{0.05};
      RealT Vfrz{0.7};
      RealT Rc{ZERO<RealT>};
      RealT Xc{ZERO<RealT>};
      RealT Kc{ONE<RealT>};
      RealT dbdlow{ZERO<RealT>};
      RealT dbdupper{ZERO<RealT>};
      RealT emax{ONE<RealT>};
      RealT emin{-ONE<RealT>};
      RealT Kp{10.0};
      RealT Ki{10.0};
      RealT Qmax{ONE<RealT>};
      RealT Qmin{-ONE<RealT>};
      RealT Tft{ZERO<RealT>};
      RealT Tfv{3.0};
      RealT Tp{ZERO<RealT>};
      RealT fdbd1{ZERO<RealT>};
      RealT fdbd2{ZERO<RealT>};
      RealT Ddn{20.0};
      RealT Dup{ZERO<RealT>};
      RealT femax{ONE<RealT>};
      RealT femin{-ONE<RealT>};
      RealT Kpg{10.0};
      RealT Kig{10.0};
      RealT Pmax{TWO<RealT>};
      RealT Pmin{ZERO<RealT>};
      RealT Tlag{3.0};
      RealT V{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::mva, &RepcaData::mva},
                                             Parameter{Parameters::VcompFlag, &RepcaData::VcompFlag},
                                             Parameter{Parameters::RefFlag, &RepcaData::RefFlag},
                                             Parameter{Parameters::Freqflag, &RepcaData::Freqflag},
                                             Parameter{Parameters::Tfltr, &RepcaData::Tfltr},
                                             Parameter{Parameters::Vfrz, &RepcaData::Vfrz},
                                             Parameter{Parameters::Rc, &RepcaData::Rc},
                                             Parameter{Parameters::Xc, &RepcaData::Xc},
                                             Parameter{Parameters::Kc, &RepcaData::Kc},
                                             Parameter{Parameters::dbdlow, &RepcaData::dbdlow},
                                             Parameter{Parameters::dbdupper, &RepcaData::dbdupper},
                                             Parameter{Parameters::emax, &RepcaData::emax},
                                             Parameter{Parameters::emin, &RepcaData::emin},
                                             Parameter{Parameters::Kp, &RepcaData::Kp},
                                             Parameter{Parameters::Ki, &RepcaData::Ki},
                                             Parameter{Parameters::Qmax, &RepcaData::Qmax},
                                             Parameter{Parameters::Qmin, &RepcaData::Qmin},
                                             Parameter{Parameters::Tft, &RepcaData::Tft},
                                             Parameter{Parameters::Tfv, &RepcaData::Tfv},
                                             Parameter{Parameters::Tp, &RepcaData::Tp},
                                             Parameter{Parameters::fdbd1, &RepcaData::fdbd1},
                                             Parameter{Parameters::fdbd2, &RepcaData::fdbd2},
                                             Parameter{Parameters::Ddn, &RepcaData::Ddn},
                                             Parameter{Parameters::Dup, &RepcaData::Dup},
                                             Parameter{Parameters::femax, &RepcaData::femax},
                                             Parameter{Parameters::femin, &RepcaData::femin},
                                             Parameter{Parameters::Kpg, &RepcaData::Kpg},
                                             Parameter{Parameters::Kig, &RepcaData::Kig},
                                             Parameter{Parameters::Pmax, &RepcaData::Pmax},
                                             Parameter{Parameters::Pmin, &RepcaData::Pmin},
                                             Parameter{Parameters::Tlag, &RepcaData::Tlag},
                                             Parameter{Parameters::V, &RepcaData::V}};
    };
  } // namespace EMT
} // namespace GridKit
