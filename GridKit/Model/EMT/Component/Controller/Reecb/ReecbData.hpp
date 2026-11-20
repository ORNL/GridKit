/**
 * @file ReecbData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT renewable energy electrical control.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Reecb`
    enum class ReecbInputs : size_t
    {
      vt,     ///< \f$V_t\f$ Terminal voltage, line-to-line RMS [V]
      p,      ///< \f$P\f$ Active power into the grid [W]
      q,      ///< \f$Q\f$ Reactive power into the grid [var]
      qext,   ///< \f$Q^\mathrm{ext}\f$ Reactive-power command on the system power base, or voltage command [p.u.]
      pfaref, ///< \f$\phi^\mathrm{ref}\f$ Power-factor angle reference [rad]
      pref,   ///< \f$P^\mathrm{ref}\f$ Active-power command on the system power base [p.u.]
    };

    /// Parameters of a `Reecb`
    enum class ReecbParameters : size_t
    {
      PfFlag, ///< \f$s_\mathrm{pf}\f$ Power-factor control selector
      VFlag,  ///< \f$s_V\f$ Voltage-reference selector of the Volt/VAr control
      QFlag,  ///< \f$s_Q\f$ Reactive-path selector
      Pqflag, ///< \f$s_\mathrm{pq}\f$ Current-priority selector
      Trv,    ///< \f$T_\mathrm{rv}\f$ Voltage-measurement filter time constant [s]
      Tp,     ///< \f$T_\mathrm{p}\f$ Electrical-power measurement filter time constant [s]
      Vref0,  ///< \f$V^\mathrm{ref}\f$ Voltage reference of the reactive-current injection [p.u.]
      Vdip,   ///< \f$V_\mathrm{dip}\f$ Low-voltage threshold for the voltage-band gate [p.u.]
      Vup,    ///< \f$V_\mathrm{up}\f$ High-voltage threshold for the voltage-band gate [p.u.]
      dbd1,   ///< \f$D_1^\mathrm{db}\f$ Lower deadband threshold for voltage-error response [p.u.]
      dbd2,   ///< \f$D_2^\mathrm{db}\f$ Upper deadband threshold for voltage-error response [p.u.]
      kqv,    ///< \f$K_\mathrm{qv}\f$ Reactive-current injection gain [p.u.]
      Iql1,   ///< \f$I_{q,\mathrm{inj}}^{\min}\f$ Minimum reactive-current injection [p.u.]
      Iqh1,   ///< \f$I_{q,\mathrm{inj}}^{\max}\f$ Maximum reactive-current injection [p.u.]
      Qmax,   ///< \f$Q^{\max}\f$ Maximum reactive-power control output [p.u.]
      Qmin,   ///< \f$Q^{\min}\f$ Minimum reactive-power control output [p.u.]
      Kqp,    ///< \f$K_\mathrm{qp}\f$ Reactive-power controller proportional gain [p.u.]
      Kqi,    ///< \f$K_\mathrm{qi}\f$ Reactive-power controller integral gain [p.u./s]
      Vmax,   ///< \f$V^{\max}\f$ Maximum voltage-control output [p.u.]
      Vmin,   ///< \f$V^{\min}\f$ Minimum voltage-control output [p.u.]
      Kvp,    ///< \f$K_\mathrm{vp}\f$ Voltage controller proportional gain [p.u.]
      Kvi,    ///< \f$K_\mathrm{vi}\f$ Voltage controller integral gain [p.u./s]
      Tiq,    ///< \f$T_\mathrm{iq}\f$ Reactive-current command lag time constant [s]
      Tpord,  ///< \f$T_\mathrm{pord}\f$ Active-power order filter time constant [s]
      dPmax,  ///< \f$R_P^{\max}\f$ Positive active-power order ramp-rate limit [p.u./s]
      dPmin,  ///< \f$R_P^{\min}\f$ Negative active-power order ramp-rate limit [p.u./s]
      Pmax,   ///< \f$P^{\max}\f$ Maximum active-power order [p.u.]
      Pmin,   ///< \f$P^{\min}\f$ Minimum active-power order [p.u.]
      Imax,   ///< \f$I^{\max}\f$ Maximum converter current [p.u.]
      V,      ///< \f$V\f$ Rated line-to-line RMS voltage [V]
      mva,    ///< \f$S\f$ Three-phase device rating [MVA]
    };

    /// Data of a `Reecb`
    template <typename RealT, typename IdxT>
    struct ReecbData
    {
      using Parameters = ReecbParameters;

      bool  PfFlag{false};
      bool  VFlag{false};
      bool  QFlag{false};
      bool  Pqflag{false};
      RealT Trv{0.02};
      RealT Tp{ZERO<RealT>};
      RealT Vref0{ZERO<RealT>};
      RealT Vdip{0.85};
      RealT Vup{1.15};
      RealT dbd1{ZERO<RealT>};
      RealT dbd2{ZERO<RealT>};
      RealT kqv{5.0};
      RealT Iql1{-1.1};
      RealT Iqh1{1.1};
      RealT Qmax{0.436};
      RealT Qmin{-0.436};
      RealT Kqp{ZERO<RealT>};
      RealT Kqi{0.1};
      RealT Vmax{1.1};
      RealT Vmin{0.9};
      RealT Kvp{18.0};
      RealT Kvi{5.0};
      RealT Tiq{0.02};
      RealT Tpord{0.02};
      RealT dPmax{99.0};
      RealT dPmin{-99.0};
      RealT Pmax{ONE<RealT>};
      RealT Pmin{ZERO<RealT>};
      RealT Imax{1.3};
      RealT V{ZERO<RealT>};
      RealT mva{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::PfFlag, &ReecbData::PfFlag},
                                             Parameter{Parameters::VFlag, &ReecbData::VFlag},
                                             Parameter{Parameters::QFlag, &ReecbData::QFlag},
                                             Parameter{Parameters::Pqflag, &ReecbData::Pqflag},
                                             Parameter{Parameters::Trv, &ReecbData::Trv},
                                             Parameter{Parameters::Tp, &ReecbData::Tp},
                                             Parameter{Parameters::Vref0, &ReecbData::Vref0},
                                             Parameter{Parameters::Vdip, &ReecbData::Vdip},
                                             Parameter{Parameters::Vup, &ReecbData::Vup},
                                             Parameter{Parameters::dbd1, &ReecbData::dbd1},
                                             Parameter{Parameters::dbd2, &ReecbData::dbd2},
                                             Parameter{Parameters::kqv, &ReecbData::kqv},
                                             Parameter{Parameters::Iql1, &ReecbData::Iql1},
                                             Parameter{Parameters::Iqh1, &ReecbData::Iqh1},
                                             Parameter{Parameters::Qmax, &ReecbData::Qmax},
                                             Parameter{Parameters::Qmin, &ReecbData::Qmin},
                                             Parameter{Parameters::Kqp, &ReecbData::Kqp},
                                             Parameter{Parameters::Kqi, &ReecbData::Kqi},
                                             Parameter{Parameters::Vmax, &ReecbData::Vmax},
                                             Parameter{Parameters::Vmin, &ReecbData::Vmin},
                                             Parameter{Parameters::Kvp, &ReecbData::Kvp},
                                             Parameter{Parameters::Kvi, &ReecbData::Kvi},
                                             Parameter{Parameters::Tiq, &ReecbData::Tiq},
                                             Parameter{Parameters::Tpord, &ReecbData::Tpord},
                                             Parameter{Parameters::dPmax, &ReecbData::dPmax},
                                             Parameter{Parameters::dPmin, &ReecbData::dPmin},
                                             Parameter{Parameters::Pmax, &ReecbData::Pmax},
                                             Parameter{Parameters::Pmin, &ReecbData::Pmin},
                                             Parameter{Parameters::Imax, &ReecbData::Imax},
                                             Parameter{Parameters::V, &ReecbData::V},
                                             Parameter{Parameters::mva, &ReecbData::mva}};
    };
  } // namespace EMT
} // namespace GridKit
