/**
 * @file GastPtiData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT gas turbine governor.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `GastPti`, each of size 1 and in per unit
    enum class GastPtiInputs : size_t
    {
      omega, ///< \f$\omega\f$ Rotor speed
      pref,  ///< \f$P^\mathrm{ref}\f$ Load reference on the system power base
    };

    /// Parameters of a `GastPti`
    enum class GastPtiParameters : size_t
    {
      R,     ///< \f$R\f$ Permanent speed droop [p.u.]
      T1,    ///< \f$T_1\f$ Fuel-valve time constant [s]
      T2,    ///< \f$T_2\f$ Fuel-flow time constant [s]
      T3,    ///< \f$T_3\f$ Exhaust-temperature time constant [s]
      At,    ///< \f$A_T\f$ Ambient-temperature load limit [p.u.]
      Kt,    ///< \f$K_T\f$ Exhaust-temperature feedback gain [p.u.]
      Vmax,  ///< \f$V^{\max}\f$ Upper valve response limit [p.u.]
      Vmin,  ///< \f$V^{\min}\f$ Lower valve response limit [p.u.]
      Dturb, ///< \f$D^\mathrm{turb}\f$ Turbine damping coefficient [p.u.]
      Trate, ///< \f$T^\mathrm{rate}\f$ Turbine rating [MW]
      S,     ///< \f$S\f$ Rated three-phase apparent power of the machine [VA]
    };

    /// Data of a `GastPti`
    template <typename RealT, typename IdxT>
    struct GastPtiData
    {
      using Parameters = GastPtiParameters;

      RealT R{0.05};
      RealT T1{0.4};
      RealT T2{0.1};
      RealT T3{3.0};
      RealT At{ONE<RealT>};
      RealT Kt{TWO<RealT>};
      RealT Vmax{ONE<RealT>};
      RealT Vmin{ZERO<RealT>};
      RealT Dturb{ZERO<RealT>};
      RealT Trate{ZERO<RealT>};
      RealT S{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::R, &GastPtiData::R},
                                             Parameter{Parameters::T1, &GastPtiData::T1},
                                             Parameter{Parameters::T2, &GastPtiData::T2},
                                             Parameter{Parameters::T3, &GastPtiData::T3},
                                             Parameter{Parameters::At, &GastPtiData::At},
                                             Parameter{Parameters::Kt, &GastPtiData::Kt},
                                             Parameter{Parameters::Vmax, &GastPtiData::Vmax},
                                             Parameter{Parameters::Vmin, &GastPtiData::Vmin},
                                             Parameter{Parameters::Dturb, &GastPtiData::Dturb},
                                             Parameter{Parameters::Trate, &GastPtiData::Trate},
                                             Parameter{Parameters::S, &GastPtiData::S}};
    };
  } // namespace EMT
} // namespace GridKit
