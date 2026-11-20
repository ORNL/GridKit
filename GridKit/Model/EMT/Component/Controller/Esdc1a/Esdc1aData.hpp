/**
 * @file Esdc1aData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT IEEE type DC1A excitation system.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of an `Esdc1a`, each of size 1
    enum class Esdc1aInputs : size_t
    {
      vt,    ///< \f$V_t\f$ Terminal voltage, line-to-line RMS [V]
      omega, ///< \f$\omega\f$ Rotor speed [p.u.]
      vref,  ///< \f$V_{\mathrm{ref}}\f$ Voltage reference [p.u.]
      vs,    ///< \f$V_S\f$ Stabilizer signal [p.u.]
      vuel,  ///< \f$V_{\mathrm{uel}}\f$ Underexcitation limiter signal [p.u.]
    };

    /// Parameters of an `Esdc1a`
    enum class Esdc1aParameters : size_t
    {
      V,      ///< \f$V\f$ Rated line-to-line RMS voltage [V]
      Tr,     ///< \f$T_R\f$ Voltage transducer time constant [s]
      Ka,     ///< \f$K_A\f$ Voltage-regulator gain [p.u.]
      Ta,     ///< \f$T_A\f$ Voltage-regulator time constant [s]
      Tb,     ///< \f$T_B\f$ Input lead-lag denominator time constant [s]
      Tc,     ///< \f$T_C\f$ Input lead-lag numerator time constant [s]
      Vrmax,  ///< \f$V_R^{\max}\f$ Maximum voltage-regulator output [p.u.]
      Vrmin,  ///< \f$V_R^{\min}\f$ Minimum voltage-regulator output [p.u.]
      Ke,     ///< \f$K_E\f$ Exciter field resistance line slope margin [p.u.]
      Te,     ///< \f$T_E\f$ Exciter time constant [s]
      Kf,     ///< \f$K_F\f$ Stabilizing feedback gain [p.u.]
      Tf1,    ///< \f$T_{F1}\f$ Stabilizing feedback time constant [s]
      Spdmlt, ///< \f$s_{\mathrm{spd}}\f$ Field-voltage speed-multiplier flag
      E1,     ///< \f$E_1\f$ First saturation voltage point [p.u.]
      Se1,    ///< \f$S_E(E_1)\f$ Saturation coefficient at $E_1$ [p.u.]
      E2,     ///< \f$E_2\f$ Second saturation voltage point [p.u.]
      Se2,    ///< \f$S_E(E_2)\f$ Saturation coefficient at $E_2$ [p.u.]
      UEL,    ///< \f$I_{\mathrm{uel}}\f$ Under-excitation limiter input-routing selector
      exclim, ///< \f$s_{\mathrm{lim}}\f$ Exciter field-voltage-state lower-limit flag
    };

    /// Data of an `Esdc1a`
    template <typename RealT, typename IdxT>
    struct Esdc1aData
    {
      using Parameters = Esdc1aParameters;

      RealT V{ZERO<RealT>};
      RealT Tr{ZERO<RealT>};
      RealT Ka{40.0};
      RealT Ta{0.1};
      RealT Tb{ZERO<RealT>};
      RealT Tc{ZERO<RealT>};
      RealT Vrmax{ONE<RealT>};
      RealT Vrmin{-ONE<RealT>};
      RealT Ke{0.1};
      RealT Te{0.5};
      RealT Kf{0.05};
      RealT Tf1{0.7};
      bool  Spdmlt{false};
      RealT E1{2.8};
      RealT Se1{0.08};
      RealT E2{3.7};
      RealT Se2{0.33};
      IdxT  UEL{0};
      bool  exclim{true};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::V, &Esdc1aData::V},
                                             Parameter{Parameters::Tr, &Esdc1aData::Tr},
                                             Parameter{Parameters::Ka, &Esdc1aData::Ka},
                                             Parameter{Parameters::Ta, &Esdc1aData::Ta},
                                             Parameter{Parameters::Tb, &Esdc1aData::Tb},
                                             Parameter{Parameters::Tc, &Esdc1aData::Tc},
                                             Parameter{Parameters::Vrmax, &Esdc1aData::Vrmax},
                                             Parameter{Parameters::Vrmin, &Esdc1aData::Vrmin},
                                             Parameter{Parameters::Ke, &Esdc1aData::Ke},
                                             Parameter{Parameters::Te, &Esdc1aData::Te},
                                             Parameter{Parameters::Kf, &Esdc1aData::Kf},
                                             Parameter{Parameters::Tf1, &Esdc1aData::Tf1},
                                             Parameter{Parameters::Spdmlt, &Esdc1aData::Spdmlt},
                                             Parameter{Parameters::E1, &Esdc1aData::E1},
                                             Parameter{Parameters::Se1, &Esdc1aData::Se1},
                                             Parameter{Parameters::E2, &Esdc1aData::E2},
                                             Parameter{Parameters::Se2, &Esdc1aData::Se2},
                                             Parameter{Parameters::UEL, &Esdc1aData::UEL},
                                             Parameter{Parameters::exclim, &Esdc1aData::exclim}};
    };
  } // namespace EMT
} // namespace GridKit
