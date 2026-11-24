/**
 * @file Esst4bData.hpp
 * @brief Modeling data for the IEEE ST4B static excitation system.
 */

#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class Esst4bInputs : size_t
    {
      vt,   ///< Terminal line-to-line RMS voltage [V]
      v,    ///< Instantaneous phase-to-neutral voltages a, b, c [V]
      i,    ///< Instantaneous terminal currents injected into the grid a, b, c [A]
      ifd,  ///< Field current on the exciter base [p.u.]
      vref, ///< Voltage reference [p.u.]
      vs,   ///< Stabilizer signal [p.u.]
      vuel, ///< Underexcitation limiter signal [p.u.]
      voel, ///< Overexcitation limiter signal [p.u.]
    };

    enum class Esst4bParameters : size_t
    {
      S, ///< Machine rating [VA]
      V, ///< Rated line-to-line RMS voltage [V]
      Tr,
      Kpr,
      Kir,
      Vrmax,
      Vrmin,
      Ta,
      Kpm,
      Kim,
      Vmmax,
      Vmmin,
      Kg,
      Kp,
      Ki,
      Vbmax,
      Kc,
      Xl,
      ThetaP, ///< Compensating voltage gain angle [degrees]
      Vgmax,
    };

    template <typename RealT, typename IdxT>
    struct Esst4bData
    {
      using Parameters = Esst4bParameters;

      RealT S{};
      RealT V{};
      RealT Tr{};
      RealT Kpr{};
      RealT Kir{};
      RealT Vrmax{};
      RealT Vrmin{};
      RealT Ta{};
      RealT Kpm{};
      RealT Kim{};
      RealT Vmmax{};
      RealT Vmmin{};
      RealT Kg{};
      RealT Kp{};
      RealT Ki{};
      RealT Vbmax{};
      RealT Kc{};
      RealT Xl{};
      RealT ThetaP{};
      RealT Vgmax{};

      static constexpr std::tuple PARAMETERS{
          Parameter{Parameters::S, &Esst4bData::S},
          Parameter{Parameters::V, &Esst4bData::V},
          Parameter{Parameters::Tr, &Esst4bData::Tr},
          Parameter{Parameters::Kpr, &Esst4bData::Kpr},
          Parameter{Parameters::Kir, &Esst4bData::Kir},
          Parameter{Parameters::Vrmax, &Esst4bData::Vrmax},
          Parameter{Parameters::Vrmin, &Esst4bData::Vrmin},
          Parameter{Parameters::Ta, &Esst4bData::Ta},
          Parameter{Parameters::Kpm, &Esst4bData::Kpm},
          Parameter{Parameters::Kim, &Esst4bData::Kim},
          Parameter{Parameters::Vmmax, &Esst4bData::Vmmax},
          Parameter{Parameters::Vmmin, &Esst4bData::Vmmin},
          Parameter{Parameters::Kg, &Esst4bData::Kg},
          Parameter{Parameters::Kp, &Esst4bData::Kp},
          Parameter{Parameters::Ki, &Esst4bData::Ki},
          Parameter{Parameters::Vbmax, &Esst4bData::Vbmax},
          Parameter{Parameters::Kc, &Esst4bData::Kc},
          Parameter{Parameters::Xl, &Esst4bData::Xl},
          Parameter{Parameters::ThetaP, &Esst4bData::ThetaP},
          Parameter{Parameters::Vgmax, &Esst4bData::Vgmax}};
    };
  } // namespace EMT
} // namespace GridKit
