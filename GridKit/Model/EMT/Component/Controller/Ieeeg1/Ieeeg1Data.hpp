/**
 * @file Ieeeg1Data.hpp
 * @brief Modeling data for the IEEE type 1 steam turbine governor.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Ieeeg1Inputs : size_t
    {
      omega, ///< Rotor speed [p.u.]
      pref,  ///< Power reference [p.u.]
      paux,  ///< Auxiliary command [p.u.]
    };

    enum class Ieeeg1Parameters : size_t
    {
      S,     ///< Machine rating [VA]
      Trate, ///< Turbine rating [MW]. Zero selects the machine rating
      K,
      T1,
      T2,
      T3,
      Uo,
      Uc,
      Pmax,
      Pmin,
      T4,
      K1,
      K2,
      T5,
      K3,
      K4,
      T6,
      K5,
      K6,
      T7,
      K7,
      K8,
      db1,
      db2,
      Eps,
      Gv,  ///< Valve stroke breakpoints [p.u.]
      Pgv, ///< Steam flow at each valve stroke [p.u.]
    };

    template <typename RealT, typename IdxT>
    struct Ieeeg1Data
    {
      using Parameters = Ieeeg1Parameters;

      RealT                S{};
      RealT                Trate{};
      RealT                K{};
      RealT                T1{};
      RealT                T2{};
      RealT                T3{};
      RealT                Uo{};
      RealT                Uc{};
      RealT                Pmax{};
      RealT                Pmin{};
      RealT                T4{};
      RealT                K1{};
      RealT                K2{};
      RealT                T5{};
      RealT                K3{};
      RealT                K4{};
      RealT                T6{};
      RealT                K5{};
      RealT                K6{};
      RealT                T7{};
      RealT                K7{};
      RealT                K8{};
      RealT                db1{};
      RealT                db2{};
      RealT                Eps{};
      std::array<RealT, 6> Gv{};
      std::array<RealT, 6> Pgv{};

      static constexpr std::tuple PARAMETERS{
          Parameter{Parameters::S, &Ieeeg1Data::S},
          Parameter{Parameters::Trate, &Ieeeg1Data::Trate},
          Parameter{Parameters::K, &Ieeeg1Data::K},
          Parameter{Parameters::T1, &Ieeeg1Data::T1},
          Parameter{Parameters::T2, &Ieeeg1Data::T2},
          Parameter{Parameters::T3, &Ieeeg1Data::T3},
          Parameter{Parameters::Uo, &Ieeeg1Data::Uo},
          Parameter{Parameters::Uc, &Ieeeg1Data::Uc},
          Parameter{Parameters::Pmax, &Ieeeg1Data::Pmax},
          Parameter{Parameters::Pmin, &Ieeeg1Data::Pmin},
          Parameter{Parameters::T4, &Ieeeg1Data::T4},
          Parameter{Parameters::K1, &Ieeeg1Data::K1},
          Parameter{Parameters::K2, &Ieeeg1Data::K2},
          Parameter{Parameters::T5, &Ieeeg1Data::T5},
          Parameter{Parameters::K3, &Ieeeg1Data::K3},
          Parameter{Parameters::K4, &Ieeeg1Data::K4},
          Parameter{Parameters::T6, &Ieeeg1Data::T6},
          Parameter{Parameters::K5, &Ieeeg1Data::K5},
          Parameter{Parameters::K6, &Ieeeg1Data::K6},
          Parameter{Parameters::T7, &Ieeeg1Data::T7},
          Parameter{Parameters::K7, &Ieeeg1Data::K7},
          Parameter{Parameters::K8, &Ieeeg1Data::K8},
          Parameter{Parameters::db1, &Ieeeg1Data::db1},
          Parameter{Parameters::db2, &Ieeeg1Data::db2},
          Parameter{Parameters::Eps, &Ieeeg1Data::Eps},
          Parameter{Parameters::Gv, &Ieeeg1Data::Gv},
          Parameter{Parameters::Pgv, &Ieeeg1Data::Pgv}};
    };
  } // namespace EMT
} // namespace GridKit
