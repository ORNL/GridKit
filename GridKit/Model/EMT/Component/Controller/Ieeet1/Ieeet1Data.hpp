/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT IEEE type 1 excitation system.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of an `Ieeet1`, each of size 1
    enum class Ieeet1Inputs : size_t
    {
      vt,    ///< \f$V_t\f$ Terminal voltage, line-to-line RMS [V]
      omega, ///< \f$\omega\f$ Rotor speed [p.u.]
      vref,  ///< \f$V_{\mathrm{ref}}\f$ Voltage reference [p.u.]
      vs,    ///< \f$V_S\f$ Stabilizer signal [p.u.]
      vuel,  ///< \f$V_{\mathrm{uel}}\f$ Underexcitation limiter signal [p.u.]
      voel,  ///< \f$V_{\mathrm{oel}}\f$ Overexcitation limiter signal [p.u.]
    };

    /// Parameters of an `Ieeet1`
    enum class Ieeet1Parameters : size_t
    {
      V,       ///< \f$V\f$ Rated line-to-line RMS voltage [V]
      Tr,      ///< \f$T_R\f$ Time constant for voltage sensing [s]
      Ka,      ///< \f$K_A\f$ Coefficient for voltage regulation [p.u.]
      Ta,      ///< \f$T_A\f$ Time constant for voltage regulation [s]
      Ke,      ///< \f$K_E\f$ Exciter field resistance line slope margin [p.u.]
      Te,      ///< \f$T_E\f$ Time constant for excitation system [s]
      Kf,      ///< \f$K_F\f$ Coefficient for feedback [p.u.]
      Tf,      ///< \f$T_F\f$ Time constant for feedback [s]
      Vrmin,   ///< \f$V_R^{\min}\f$ Lower limit to voltage regulation [p.u.]
      Vrmax,   ///< \f$V_R^{\max}\f$ Upper limit to voltage regulation [p.u.]
      E1,      ///< \f$E_1\f$ First saturation voltage point [p.u.]
      E2,      ///< \f$E_2\f$ Second saturation voltage point [p.u.]
      Se1,     ///< \f$S_E(E_1)\f$ Saturation coefficient at \f$E_1\f$ [p.u.]
      Se2,     ///< \f$S_E(E_2)\f$ Saturation coefficient at \f$E_2\f$ [p.u.]
      Ispdlim, ///< \f$I_{\mathrm{spdlim}}\f$ Speed multiplier flag, 0 or 1
    };

    /// Data of an `Ieeet1`
    template <typename RealT, typename IdxT>
    struct Ieeet1Data
    {
      using Parameters = Ieeet1Parameters;

      RealT V{ZERO<RealT>};
      RealT Tr{ZERO<RealT>};
      RealT Ka{50.0};
      RealT Ta{0.04};
      RealT Ke{-0.06};
      RealT Te{0.6};
      RealT Kf{0.09};
      RealT Tf{1.46};
      RealT Vrmin{-ONE<RealT>};
      RealT Vrmax{ONE<RealT>};
      RealT E1{2.8};
      RealT E2{3.73};
      RealT Se1{0.04};
      RealT Se2{0.33};
      RealT Ispdlim{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::V, &Ieeet1Data::V},
                                             Parameter{Parameters::Tr, &Ieeet1Data::Tr},
                                             Parameter{Parameters::Ka, &Ieeet1Data::Ka},
                                             Parameter{Parameters::Ta, &Ieeet1Data::Ta},
                                             Parameter{Parameters::Ke, &Ieeet1Data::Ke},
                                             Parameter{Parameters::Te, &Ieeet1Data::Te},
                                             Parameter{Parameters::Kf, &Ieeet1Data::Kf},
                                             Parameter{Parameters::Tf, &Ieeet1Data::Tf},
                                             Parameter{Parameters::Vrmin, &Ieeet1Data::Vrmin},
                                             Parameter{Parameters::Vrmax, &Ieeet1Data::Vrmax},
                                             Parameter{Parameters::E1, &Ieeet1Data::E1},
                                             Parameter{Parameters::E2, &Ieeet1Data::E2},
                                             Parameter{Parameters::Se1, &Ieeet1Data::Se1},
                                             Parameter{Parameters::Se2, &Ieeet1Data::Se2},
                                             Parameter{Parameters::Ispdlim, &Ieeet1Data::Ispdlim}};
    };
  } // namespace EMT
} // namespace GridKit
