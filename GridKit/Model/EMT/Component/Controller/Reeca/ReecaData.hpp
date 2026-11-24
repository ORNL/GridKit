/**
 * @file ReecaData.hpp
 * @brief Modeling data for the WECC type A renewable electrical controller.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ReecaInputs : size_t
    {
      vt,     ///< Terminal line-to-line RMS voltage [V]
      p,      ///< Terminal active power injected into the grid [W]
      q,      ///< Terminal reactive power injected into the grid [var]
      qext,   ///< Reactive power reference on the system base, or direct voltage reference [p.u.]
      pfaref, ///< Power-factor angle [rad]
      pref,   ///< Power reference [p.u.]
      omega,  ///< Rotor speed [p.u.]
    };

    enum class ReecaParameters : size_t
    {
      V,      ///< Rated line-to-line RMS voltage [V]
      mva,    ///< \f$S\f$ Three-phase device rating [MVA]
      PfFlag, ///< \f$s_\mathrm{pf}\f$ Power-factor control selector
      VFlag,  ///< \f$s_V\f$ Voltage-reference selector of the Volt/VAr control
      QFlag,  ///< \f$s_Q\f$ Reactive-path selector
      Pqflag, ///< \f$s_\mathrm{pq}\f$ Current-priority selector
      Pflag,  ///< Multiply active power reference by rotor speed
      Trv,    ///< \f$T_\mathrm{rv}\f$ Voltage-measurement filter time constant [s]
      Tp,     ///< \f$T_\mathrm{p}\f$ Electrical-power measurement filter time constant [s]
      Vref0,  ///< \f$V^\mathrm{ref}\f$ Voltage reference of the reactive-current injection [p.u.]
      Vref1,  ///< Direct voltage reference bias [p.u.]
      Vdip,   ///< \f$V_\mathrm{dip}\f$ Low-voltage threshold for the voltage-band gate [p.u.]
      Vup,    ///< \f$V_\mathrm{up}\f$ High-voltage threshold for the voltage-band gate [p.u.]
      dbd1,   ///< \f$D_1^\mathrm{db}\f$ Lower deadband threshold for voltage-error response [p.u.]
      dbd2,   ///< \f$D_2^\mathrm{db}\f$ Upper deadband threshold for voltage-error response [p.u.]
      kqv,    ///< \f$K_\mathrm{qv}\f$ Reactive-current injection gain [p.u.]
      Iql1,   ///< \f$I_{q,\mathrm{inj}}^{\min}\f$ Minimum reactive-current injection [p.u.]
      Iqh1,   ///< \f$I_{q,\mathrm{inj}}^{\max}\f$ Maximum reactive-current injection [p.u.]
      Iqfrz,  ///< Reactive injection held for positive Thld [p.u.]
      Thld,   ///< Signed post-dip reactive injection hold time [s]
      Thld2,  ///< Post-dip hold time of the active current limit [s]
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
      Vq,     ///< Reactive current table voltage breakpoints [p.u.]
      Iq,     ///< Reactive current table limits [p.u.]
      Vp,     ///< Active current table voltage breakpoints [p.u.]
      Ip,     ///< Active current table limits [p.u.]
    };

    template <typename RealT, typename IdxT>
    struct ReecaData
    {
      using Parameters = ReecaParameters;

      RealT                V{ZERO<RealT>};
      RealT                mva{ZERO<RealT>};
      bool                 PfFlag{};
      bool                 VFlag{};
      bool                 QFlag{};
      bool                 Pqflag{};
      bool                 Pflag{};
      RealT                Trv{0.02};
      RealT                Tp{ZERO<RealT>};
      RealT                Vref0{ZERO<RealT>};
      RealT                Vref1{};
      RealT                Vdip{0.85};
      RealT                Vup{1.15};
      RealT                dbd1{ZERO<RealT>};
      RealT                dbd2{ZERO<RealT>};
      RealT                kqv{5.0};
      RealT                Iql1{-1.1};
      RealT                Iqh1{1.1};
      RealT                Iqfrz{};
      RealT                Thld{};
      RealT                Thld2{};
      RealT                Qmax{0.436};
      RealT                Qmin{-0.436};
      RealT                Kqp{ZERO<RealT>};
      RealT                Kqi{0.1};
      RealT                Vmax{1.1};
      RealT                Vmin{0.9};
      RealT                Kvp{18.0};
      RealT                Kvi{5.0};
      RealT                Tiq{0.02};
      RealT                Tpord{0.02};
      RealT                dPmax{99.0};
      RealT                dPmin{-99.0};
      RealT                Pmax{ONE<RealT>};
      RealT                Pmin{ZERO<RealT>};
      RealT                Imax{1.3};
      std::array<RealT, 4> Vq{};
      std::array<RealT, 4> Iq{};
      std::array<RealT, 4> Vp{};
      std::array<RealT, 4> Ip{};

      static constexpr std::tuple PARAMETERS{
          Parameter{Parameters::V, &ReecaData::V},
          Parameter{Parameters::mva, &ReecaData::mva},
          Parameter{Parameters::PfFlag, &ReecaData::PfFlag},
          Parameter{Parameters::VFlag, &ReecaData::VFlag},
          Parameter{Parameters::QFlag, &ReecaData::QFlag},
          Parameter{Parameters::Pqflag, &ReecaData::Pqflag},
          Parameter{Parameters::Pflag, &ReecaData::Pflag},
          Parameter{Parameters::Trv, &ReecaData::Trv},
          Parameter{Parameters::Tp, &ReecaData::Tp},
          Parameter{Parameters::Vref0, &ReecaData::Vref0},
          Parameter{Parameters::Vref1, &ReecaData::Vref1},
          Parameter{Parameters::Vdip, &ReecaData::Vdip},
          Parameter{Parameters::Vup, &ReecaData::Vup},
          Parameter{Parameters::dbd1, &ReecaData::dbd1},
          Parameter{Parameters::dbd2, &ReecaData::dbd2},
          Parameter{Parameters::kqv, &ReecaData::kqv},
          Parameter{Parameters::Iql1, &ReecaData::Iql1},
          Parameter{Parameters::Iqh1, &ReecaData::Iqh1},
          Parameter{Parameters::Iqfrz, &ReecaData::Iqfrz},
          Parameter{Parameters::Thld, &ReecaData::Thld},
          Parameter{Parameters::Thld2, &ReecaData::Thld2},
          Parameter{Parameters::Qmax, &ReecaData::Qmax},
          Parameter{Parameters::Qmin, &ReecaData::Qmin},
          Parameter{Parameters::Kqp, &ReecaData::Kqp},
          Parameter{Parameters::Kqi, &ReecaData::Kqi},
          Parameter{Parameters::Vmax, &ReecaData::Vmax},
          Parameter{Parameters::Vmin, &ReecaData::Vmin},
          Parameter{Parameters::Kvp, &ReecaData::Kvp},
          Parameter{Parameters::Kvi, &ReecaData::Kvi},
          Parameter{Parameters::Tiq, &ReecaData::Tiq},
          Parameter{Parameters::Tpord, &ReecaData::Tpord},
          Parameter{Parameters::dPmax, &ReecaData::dPmax},
          Parameter{Parameters::dPmin, &ReecaData::dPmin},
          Parameter{Parameters::Pmax, &ReecaData::Pmax},
          Parameter{Parameters::Pmin, &ReecaData::Pmin},
          Parameter{Parameters::Imax, &ReecaData::Imax},
          Parameter{Parameters::Vq, &ReecaData::Vq},
          Parameter{Parameters::Iq, &ReecaData::Iq},
          Parameter{Parameters::Vp, &ReecaData::Vp},
          Parameter{Parameters::Ip, &ReecaData::Ip}};
    };
  } // namespace EMT
} // namespace GridKit
