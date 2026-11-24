/**
 * @file TransformerData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT transformers.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Terminals of a `Transformer`, each of size 3
    enum class TransformerTerminals : size_t
    {
      bus1, ///< Terminal 1
      bus2, ///< Terminal 2
    };

    /// Parameters of a `Transformer`
    enum class TransformerParameters : size_t
    {
      N,     ///< \f$N\f$ Number of phases
      S,     ///< \f$S_\mathrm{b}\f$ Rated three-phase apparent power [VA]
      V1,    ///< \f$V_1\f$ Rated line-to-line RMS voltage of winding 1 [V]
      V2,    ///< \f$V_2\f$ Rated line-to-line RMS voltage of winding 2 [V]
      f,     ///< \f$f_\mathrm{b}\f$ Rated frequency [Hz]
      P1,    ///< \f$\mathbf{P}_1\f$ Terminal 1 connection map, N by N
      P2,    ///< \f$\mathbf{P}_2\f$ Terminal 2 connection map, N by N
      tap,   ///< \f$\tau\f$ Off-nominal ratio on winding 1 [p.u.]
      phase, ///< \f$\theta\f$ Phase shift of winding 1 [rad]
      R,     ///< \f$R\f$ Short-circuit resistance [p.u.]
      X,     ///< \f$X\f$ Short-circuit reactance [p.u.]
      I0,    ///< \f$I_0\f$ No-load current at rated voltage [p.u.]
      P0,    ///< \f$P_0\f$ No-load loss at rated voltage [W]
      knee,  ///< \f$\psi_\mathrm{K}\f$ Knee flux linkage [p.u.]
      Lsat,  ///< \f$L_\mathrm{sat}\f$ Terminal saturation inductance [p.u.]
      split, ///< \f$\beta\f$ Winding 1 share of the magnetizing branch
    };

    /// Data of a `Transformer`
    template <typename RealT, typename IdxT>
    struct TransformerData
    {
      using Parameters = TransformerParameters;

      IdxT          N{3};
      RealT         S{ZERO<RealT>};
      RealT         V1{ZERO<RealT>};
      RealT         V2{ZERO<RealT>};
      RealT         f{ZERO<RealT>};
      Matrix<RealT> P1;
      Matrix<RealT> P2;
      RealT         tap{ONE<RealT>};
      RealT         phase{ZERO<RealT>};
      RealT         R{ZERO<RealT>};
      RealT         X{ZERO<RealT>};
      RealT         I0{ZERO<RealT>};
      RealT         P0{ZERO<RealT>};
      RealT         knee{1.2};
      RealT         Lsat{ZERO<RealT>}; ///< Zero selects 2 X
      RealT         split{HALF<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &TransformerData::N},
                                             Parameter{Parameters::S, &TransformerData::S},
                                             Parameter{Parameters::V1, &TransformerData::V1},
                                             Parameter{Parameters::V2, &TransformerData::V2},
                                             Parameter{Parameters::f, &TransformerData::f},
                                             Parameter{Parameters::P1, &TransformerData::P1},
                                             Parameter{Parameters::P2, &TransformerData::P2},
                                             Parameter{Parameters::tap, &TransformerData::tap},
                                             Parameter{Parameters::phase, &TransformerData::phase},
                                             Parameter{Parameters::R, &TransformerData::R},
                                             Parameter{Parameters::X, &TransformerData::X},
                                             Parameter{Parameters::I0, &TransformerData::I0},
                                             Parameter{Parameters::P0, &TransformerData::P0},
                                             Parameter{Parameters::knee, &TransformerData::knee},
                                             Parameter{Parameters::Lsat, &TransformerData::Lsat},
                                             Parameter{Parameters::split, &TransformerData::split}};
    };
  } // namespace EMT
} // namespace GridKit
