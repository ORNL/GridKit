/**
 * @file Ggov1Data.hpp
 * @brief Modeling data for the general turbine governor.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Ggov1Inputs : size_t
    {
      omega,  ///< Rotor speed [p.u.]
      p,      ///< Terminal active power injected into the grid [W]
      pref,   ///< Power reference [p.u.]
      paux,   ///< Auxiliary command [p.u.]
      pmwset, ///< Supervisory electrical power reference on the turbine rating [p.u.]
    };

    enum class Ggov1Parameters : size_t
    {
      S,       ///< Machine rating [VA]
      Trate,   ///< Turbine rating [MW]. Zero selects the machine rating
      Rselect, ///< Droop feedback selector
      Flag,    ///< Speed-dependent fuel flow
      R,
      Tpelec,
      Maxerr,
      Minerr,
      Kpgov,
      Kigov,
      Kdgov,
      Tdgov,
      Vmax,
      Vmin,
      Tact,
      Kturb,
      Wfnl,
      Tb,
      Tc,
      Teng,
      Tfload,
      Kpload,
      Kiload,
      Ldref,
      Dm,
      Ropen,
      Rclose,
      Kimw,
      Aset,
      Ka,
      Ta,
      Db,
      Tsa,
      Tsb,
      Tfsr,    ///< Acceleration limiter fuel-command memory [s]
      Fsrtmax, ///< Upper load-controller fuel limit [p.u.]
    };

    template <typename RealT, typename IdxT>
    struct Ggov1Data
    {
      using Parameters = Ggov1Parameters;

      RealT S{};
      RealT Trate{};
      int   Rselect{1};
      bool  Flag{};
      RealT R{};
      RealT Tpelec{};
      RealT Maxerr{};
      RealT Minerr{};
      RealT Kpgov{};
      RealT Kigov{};
      RealT Kdgov{};
      RealT Tdgov{};
      RealT Vmax{};
      RealT Vmin{};
      RealT Tact{};
      RealT Kturb{};
      RealT Wfnl{};
      RealT Tb{};
      RealT Tc{};
      RealT Teng{};
      RealT Tfload{};
      RealT Kpload{};
      RealT Kiload{};
      RealT Ldref{};
      RealT Dm{};
      RealT Ropen{};
      RealT Rclose{};
      RealT Kimw{};
      RealT Aset{};
      RealT Ka{};
      RealT Ta{};
      RealT Db{};
      RealT Tsa{};
      RealT Tsb{};
      RealT Tfsr{1e-3};
      RealT Fsrtmax{1};

      static constexpr std::tuple PARAMETERS{
          Parameter{Parameters::S, &Ggov1Data::S},
          Parameter{Parameters::Trate, &Ggov1Data::Trate},
          Parameter{Parameters::Rselect, &Ggov1Data::Rselect},
          Parameter{Parameters::Flag, &Ggov1Data::Flag},
          Parameter{Parameters::R, &Ggov1Data::R},
          Parameter{Parameters::Tpelec, &Ggov1Data::Tpelec},
          Parameter{Parameters::Maxerr, &Ggov1Data::Maxerr},
          Parameter{Parameters::Minerr, &Ggov1Data::Minerr},
          Parameter{Parameters::Kpgov, &Ggov1Data::Kpgov},
          Parameter{Parameters::Kigov, &Ggov1Data::Kigov},
          Parameter{Parameters::Kdgov, &Ggov1Data::Kdgov},
          Parameter{Parameters::Tdgov, &Ggov1Data::Tdgov},
          Parameter{Parameters::Vmax, &Ggov1Data::Vmax},
          Parameter{Parameters::Vmin, &Ggov1Data::Vmin},
          Parameter{Parameters::Tact, &Ggov1Data::Tact},
          Parameter{Parameters::Kturb, &Ggov1Data::Kturb},
          Parameter{Parameters::Wfnl, &Ggov1Data::Wfnl},
          Parameter{Parameters::Tb, &Ggov1Data::Tb},
          Parameter{Parameters::Tc, &Ggov1Data::Tc},
          Parameter{Parameters::Teng, &Ggov1Data::Teng},
          Parameter{Parameters::Tfload, &Ggov1Data::Tfload},
          Parameter{Parameters::Kpload, &Ggov1Data::Kpload},
          Parameter{Parameters::Kiload, &Ggov1Data::Kiload},
          Parameter{Parameters::Ldref, &Ggov1Data::Ldref},
          Parameter{Parameters::Dm, &Ggov1Data::Dm},
          Parameter{Parameters::Ropen, &Ggov1Data::Ropen},
          Parameter{Parameters::Rclose, &Ggov1Data::Rclose},
          Parameter{Parameters::Kimw, &Ggov1Data::Kimw},
          Parameter{Parameters::Aset, &Ggov1Data::Aset},
          Parameter{Parameters::Ka, &Ggov1Data::Ka},
          Parameter{Parameters::Ta, &Ggov1Data::Ta},
          Parameter{Parameters::Db, &Ggov1Data::Db},
          Parameter{Parameters::Tsa, &Ggov1Data::Tsa},
          Parameter{Parameters::Tsb, &Ggov1Data::Tsb},
          Parameter{Parameters::Tfsr, &Ggov1Data::Tfsr},
          Parameter{Parameters::Fsrtmax, &Ggov1Data::Fsrtmax}};
    };
  } // namespace EMT
} // namespace GridKit
