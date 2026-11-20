/**
 * @file Tgov1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT steam turbine governor.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Tgov1`, each of size 1 and in per unit
    enum class Tgov1Inputs : size_t
    {
      omega, ///< \f$\omega\f$ Rotor speed
      pref,  ///< \f$P_\mathrm{ref}\f$ Load reference on the turbine rating
    };

    /// Parameters of a `Tgov1`
    enum class Tgov1Parameters : size_t
    {
      Trate, ///< \f$T_{\mathrm{rate}}\f$ Turbine rating [MW]
      R,     ///< \f$R\f$ Permanent droop [p.u.]
      Pvmin, ///< \f$P_\mathrm{v}^{\min}\f$ Minimum valve position [p.u.]
      Pvmax, ///< \f$P_\mathrm{v}^{\max}\f$ Maximum valve position [p.u.]
      T1,    ///< \f$T_1\f$ Steam-bowl time constant [s]
      T2,    ///< \f$T_2\f$ Turbine numerator time constant [s]
      T3,    ///< \f$T_3\f$ Reheater time constant [s]
      Dt,    ///< \f$D_\mathrm{t}\f$ Turbine damping coefficient [p.u.]
      S,     ///< \f$S\f$ Rated three-phase apparent power of the machine [VA]
    };

    /// Data of a `Tgov1`
    template <typename RealT, typename IdxT>
    struct Tgov1Data
    {
      using Parameters = Tgov1Parameters;

      RealT Trate{100.0};
      RealT R{0.05};
      RealT Pvmin{ZERO<RealT>};
      RealT Pvmax{ONE<RealT>};
      RealT T1{0.5};
      RealT T2{2.5};
      RealT T3{7.5};
      RealT Dt{ZERO<RealT>};
      RealT S{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::Trate, &Tgov1Data::Trate},
                                             Parameter{Parameters::R, &Tgov1Data::R},
                                             Parameter{Parameters::Pvmin, &Tgov1Data::Pvmin},
                                             Parameter{Parameters::Pvmax, &Tgov1Data::Pvmax},
                                             Parameter{Parameters::T1, &Tgov1Data::T1},
                                             Parameter{Parameters::T2, &Tgov1Data::T2},
                                             Parameter{Parameters::T3, &Tgov1Data::T3},
                                             Parameter{Parameters::Dt, &Tgov1Data::Dt},
                                             Parameter{Parameters::S, &Tgov1Data::S}};
    };
  } // namespace EMT
} // namespace GridKit
