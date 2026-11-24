/**
 * @file TransformerThreeWindingData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT three-winding transformers.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Terminals of a `TransformerThreeWinding`, each of size 3
    enum class TransformerThreeWindingTerminals : size_t
    {
      bus1, ///< Terminal 1
      bus2, ///< Terminal 2
      bus3, ///< Terminal 3
    };

    /// Parameters of a `TransformerThreeWinding`
    enum class TransformerThreeWindingParameters : size_t
    {
      N,     ///< \f$N\f$ Number of phases
      S,     ///< \f$S_\mathrm{b}\f$ Rated three-phase apparent power [VA]
      V1,    ///< \f$V_1\f$ Rated line-to-line RMS voltage of winding 1 [V]
      V2,    ///< \f$V_2\f$ Rated line-to-line RMS voltage of winding 2 [V]
      V3,    ///< \f$V_3\f$ Rated line-to-line RMS voltage of winding 3 [V]
      f,     ///< \f$f_\mathrm{b}\f$ Rated frequency [Hz]
      P1,    ///< \f$\mathbf{P}_1\f$ Terminal 1 connection map, N by N
      P2,    ///< \f$\mathbf{P}_2\f$ Terminal 2 connection map, N by N
      P3,    ///< \f$\mathbf{P}_3\f$ Terminal 3 connection map, N by N
      tap,   ///< \f$\tau\f$ Off-nominal ratio on winding 1 [p.u.]
      phase, ///< \f$\theta\f$ Phase shift of winding 1 [rad]
      R12,   ///< \f$R_{12}\f$ Short-circuit resistance of windings 1 and 2 [p.u.]
      R23,   ///< \f$R_{23}\f$ Short-circuit resistance of windings 2 and 3 [p.u.]
      R31,   ///< \f$R_{31}\f$ Short-circuit resistance of windings 3 and 1 [p.u.]
      X12,   ///< \f$X_{12}\f$ Short-circuit reactance of windings 1 and 2 [p.u.]
      X23,   ///< \f$X_{23}\f$ Short-circuit reactance of windings 2 and 3 [p.u.]
      X31,   ///< \f$X_{31}\f$ Short-circuit reactance of windings 3 and 1 [p.u.]
      I0,    ///< \f$I_0\f$ No-load current at rated voltage [p.u.]
      P0,    ///< \f$P_0\f$ No-load loss at rated voltage [W]
      knee,  ///< \f$\psi_\mathrm{K}\f$ Knee flux linkage [p.u.]
      Lsat,  ///< \f$L_\mathrm{sat}\f$ Terminal saturation inductance [p.u.]
    };

    /// Data of a `TransformerThreeWinding`
    template <typename RealT, typename IdxT>
    struct TransformerThreeWindingData
    {
      using Parameters = TransformerThreeWindingParameters;

      IdxT          N{3};
      RealT         S{ZERO<RealT>};
      RealT         V1{ZERO<RealT>};
      RealT         V2{ZERO<RealT>};
      RealT         V3{ZERO<RealT>};
      RealT         f{ZERO<RealT>};
      Matrix<RealT> P1;
      Matrix<RealT> P2;
      Matrix<RealT> P3;
      RealT         tap{ONE<RealT>};
      RealT         phase{ZERO<RealT>};
      RealT         R12{ZERO<RealT>};
      RealT         R23{ZERO<RealT>};
      RealT         R31{ZERO<RealT>};
      RealT         X12{ZERO<RealT>};
      RealT         X23{ZERO<RealT>};
      RealT         X31{ZERO<RealT>};
      RealT         I0{ZERO<RealT>};
      RealT         P0{ZERO<RealT>};
      RealT         knee{1.2};
      RealT         Lsat{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &TransformerThreeWindingData::N},
                                             Parameter{Parameters::S, &TransformerThreeWindingData::S},
                                             Parameter{Parameters::V1, &TransformerThreeWindingData::V1},
                                             Parameter{Parameters::V2, &TransformerThreeWindingData::V2},
                                             Parameter{Parameters::V3, &TransformerThreeWindingData::V3},
                                             Parameter{Parameters::f, &TransformerThreeWindingData::f},
                                             Parameter{Parameters::P1, &TransformerThreeWindingData::P1},
                                             Parameter{Parameters::P2, &TransformerThreeWindingData::P2},
                                             Parameter{Parameters::P3, &TransformerThreeWindingData::P3},
                                             Parameter{Parameters::tap, &TransformerThreeWindingData::tap},
                                             Parameter{Parameters::phase, &TransformerThreeWindingData::phase},
                                             Parameter{Parameters::R12, &TransformerThreeWindingData::R12},
                                             Parameter{Parameters::R23, &TransformerThreeWindingData::R23},
                                             Parameter{Parameters::R31, &TransformerThreeWindingData::R31},
                                             Parameter{Parameters::X12, &TransformerThreeWindingData::X12},
                                             Parameter{Parameters::X23, &TransformerThreeWindingData::X23},
                                             Parameter{Parameters::X31, &TransformerThreeWindingData::X31},
                                             Parameter{Parameters::I0, &TransformerThreeWindingData::I0},
                                             Parameter{Parameters::P0, &TransformerThreeWindingData::P0},
                                             Parameter{Parameters::knee, &TransformerThreeWindingData::knee},
                                             Parameter{Parameters::Lsat, &TransformerThreeWindingData::Lsat}};
    };
  } // namespace EMT
} // namespace GridKit
