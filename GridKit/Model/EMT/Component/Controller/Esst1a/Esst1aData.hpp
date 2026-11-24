/**
 * @file Esst1aData.hpp
 * @brief Modeling data for the IEEE ST1A static excitation system.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Esst1aInputs : size_t
    {
      vt,   ///< Terminal line-to-line RMS voltage [V]
      ifd,  ///< Field current on the exciter base [p.u.]
      vref, ///< Voltage reference [p.u.]
      vs,   ///< Stabilizer signal [p.u.]
      vuel, ///< Underexcitation limiter signal [p.u.]
      voel, ///< Overexcitation limiter signal [p.u.]
    };

    enum class Esst1aParameters : size_t
    {
      V, ///< Rated line-to-line RMS voltage [V]
      Tr,
      Vimax,
      Vimin,
      Tc,
      Tb,
      Ka,
      Ta,
      Vrmax,
      Vrmin,
      Kc,
      Kf,
      Tf,
      UEL, ///< Underexcitation limiter input location
      VOS, ///< Stabilizer input location
      Tc1,
      Tb1,
      Vamax,
      Vamin,
      Klr,
      Ilr,
    };

    template <typename RealT, typename IdxT>
    struct Esst1aData
    {
      using Parameters = Esst1aParameters;

      RealT V{};
      RealT Tr{};
      RealT Vimax{};
      RealT Vimin{};
      RealT Tc{};
      RealT Tb{};
      RealT Ka{};
      RealT Ta{};
      RealT Vrmax{};
      RealT Vrmin{};
      RealT Kc{};
      RealT Kf{};
      RealT Tf{};
      int   UEL{1};
      int   VOS{1};
      RealT Tc1{};
      RealT Tb1{};
      RealT Vamax{};
      RealT Vamin{};
      RealT Klr{};
      RealT Ilr{};

      static constexpr std::tuple PARAMETERS{
          Parameter{Parameters::V, &Esst1aData::V},
          Parameter{Parameters::Tr, &Esst1aData::Tr},
          Parameter{Parameters::Vimax, &Esst1aData::Vimax},
          Parameter{Parameters::Vimin, &Esst1aData::Vimin},
          Parameter{Parameters::Tc, &Esst1aData::Tc},
          Parameter{Parameters::Tb, &Esst1aData::Tb},
          Parameter{Parameters::Ka, &Esst1aData::Ka},
          Parameter{Parameters::Ta, &Esst1aData::Ta},
          Parameter{Parameters::Vrmax, &Esst1aData::Vrmax},
          Parameter{Parameters::Vrmin, &Esst1aData::Vrmin},
          Parameter{Parameters::Kc, &Esst1aData::Kc},
          Parameter{Parameters::Kf, &Esst1aData::Kf},
          Parameter{Parameters::Tf, &Esst1aData::Tf},
          Parameter{Parameters::UEL, &Esst1aData::UEL},
          Parameter{Parameters::VOS, &Esst1aData::VOS},
          Parameter{Parameters::Tc1, &Esst1aData::Tc1},
          Parameter{Parameters::Tb1, &Esst1aData::Tb1},
          Parameter{Parameters::Vamax, &Esst1aData::Vamax},
          Parameter{Parameters::Vamin, &Esst1aData::Vamin},
          Parameter{Parameters::Klr, &Esst1aData::Klr},
          Parameter{Parameters::Ilr, &Esst1aData::Ilr}};
    };
  } // namespace EMT
} // namespace GridKit
