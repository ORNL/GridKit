/**
 * @file Exst1PtiData.hpp
 * @brief Modeling data for the PTI static excitation system.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Exst1PtiInputs : size_t
    {
      vt,   ///< Terminal line-to-line RMS voltage [V]
      ifd,  ///< Field current on the exciter base [p.u.]
      vref, ///< Voltage reference [p.u.]
      vs,   ///< Stabilizer signal [p.u.]
      vuel, ///< Underexcitation limiter signal [p.u.]
      voel, ///< Overexcitation limiter signal [p.u.]
    };

    enum class Exst1PtiParameters : size_t
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
    };

    template <typename RealT, typename IdxT>
    struct Exst1PtiData
    {
      using Parameters = Exst1PtiParameters;

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

      static constexpr std::tuple PARAMETERS{
          Parameter{Parameters::V, &Exst1PtiData::V},
          Parameter{Parameters::Tr, &Exst1PtiData::Tr},
          Parameter{Parameters::Vimax, &Exst1PtiData::Vimax},
          Parameter{Parameters::Vimin, &Exst1PtiData::Vimin},
          Parameter{Parameters::Tc, &Exst1PtiData::Tc},
          Parameter{Parameters::Tb, &Exst1PtiData::Tb},
          Parameter{Parameters::Ka, &Exst1PtiData::Ka},
          Parameter{Parameters::Ta, &Exst1PtiData::Ta},
          Parameter{Parameters::Vrmax, &Exst1PtiData::Vrmax},
          Parameter{Parameters::Vrmin, &Exst1PtiData::Vrmin},
          Parameter{Parameters::Kc, &Exst1PtiData::Kc},
          Parameter{Parameters::Kf, &Exst1PtiData::Kf},
          Parameter{Parameters::Tf, &Exst1PtiData::Tf}};
    };
  } // namespace EMT
} // namespace GridKit
