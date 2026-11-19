/**
 * @file IeeestData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT IEEE power system stabilizer.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of an `Ieeest`
    enum class IeeestInputs : size_t
    {
      u, ///< \f$u\f$ Stabilizer input
    };

    /// Parameters of an `Ieeest`
    enum class IeeestParameters : size_t
    {
      A1,     ///< \f$A_1\f$ Notch denominator coefficient [s]
      A2,     ///< \f$A_2\f$ Notch denominator coefficient [s^2]
      A3,     ///< \f$A_3\f$ Notch denominator coefficient [s]
      A4,     ///< \f$A_4\f$ Notch denominator coefficient [s^2]
      A5,     ///< \f$A_5\f$ Notch numerator coefficient [s]
      A6,     ///< \f$A_6\f$ Notch numerator coefficient [s^2]
      T1,     ///< \f$T_1\f$ Lead-lag 1 numerator time constant [s]
      T2,     ///< \f$T_2\f$ Lead-lag 1 denominator time constant [s]
      T3,     ///< \f$T_3\f$ Lead-lag 2 numerator time constant [s]
      T4,     ///< \f$T_4\f$ Lead-lag 2 denominator time constant [s]
      T5,     ///< \f$T_5\f$ Washout numerator time constant [s]
      T6,     ///< \f$T_6\f$ Washout denominator time constant [s]
      Ks,     ///< \f$K_\mathrm{s}\f$ Stabilizer gain [p.u.]
      Lsmin,  ///< \f$L_\mathrm{s}^{\min}\f$ Minimum stabilizer output limit [p.u.]
      Lsmax,  ///< \f$L_\mathrm{s}^{\max}\f$ Maximum stabilizer output limit [p.u.]
      Vcl,    ///< \f$V_\mathrm{cl}\f$ Lower cutout voltage, which is not supported [p.u.]
      Vcu,    ///< \f$V_\mathrm{cu}\f$ Upper cutout voltage, which is not supported [p.u.]
      Tdelay, ///< \f$T_\mathrm{delay}\f$ Input delay, which is not supported [s]
    };

    /// Data of an `Ieeest`
    template <typename RealT, typename IdxT>
    struct IeeestData
    {
      using Parameters = IeeestParameters;

      RealT A1{ZERO<RealT>};
      RealT A2{ZERO<RealT>};
      RealT A3{ZERO<RealT>};
      RealT A4{ZERO<RealT>};
      RealT A5{ZERO<RealT>};
      RealT A6{ZERO<RealT>};
      RealT T1{ZERO<RealT>};
      RealT T2{ONE<RealT>};
      RealT T3{ZERO<RealT>};
      RealT T4{ONE<RealT>};
      RealT T5{ZERO<RealT>};
      RealT T6{ONE<RealT>};
      RealT Ks{ONE<RealT>};
      RealT Lsmin{-0.1};
      RealT Lsmax{0.1};
      RealT Vcl{ZERO<RealT>};
      RealT Vcu{ZERO<RealT>};
      RealT Tdelay{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::A1, &IeeestData::A1},
                                             Parameter{Parameters::A2, &IeeestData::A2},
                                             Parameter{Parameters::A3, &IeeestData::A3},
                                             Parameter{Parameters::A4, &IeeestData::A4},
                                             Parameter{Parameters::A5, &IeeestData::A5},
                                             Parameter{Parameters::A6, &IeeestData::A6},
                                             Parameter{Parameters::T1, &IeeestData::T1},
                                             Parameter{Parameters::T2, &IeeestData::T2},
                                             Parameter{Parameters::T3, &IeeestData::T3},
                                             Parameter{Parameters::T4, &IeeestData::T4},
                                             Parameter{Parameters::T5, &IeeestData::T5},
                                             Parameter{Parameters::T6, &IeeestData::T6},
                                             Parameter{Parameters::Ks, &IeeestData::Ks},
                                             Parameter{Parameters::Lsmin, &IeeestData::Lsmin},
                                             Parameter{Parameters::Lsmax, &IeeestData::Lsmax},
                                             Parameter{Parameters::Vcl, &IeeestData::Vcl},
                                             Parameter{Parameters::Vcu, &IeeestData::Vcu},
                                             Parameter{Parameters::Tdelay, &IeeestData::Tdelay}};
    };
  } // namespace EMT
} // namespace GridKit
