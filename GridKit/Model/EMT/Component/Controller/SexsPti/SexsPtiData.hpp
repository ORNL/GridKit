/**
 * @file SexsPtiData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for an EMT simplified excitation system.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `SexsPti`, each of size 1
    enum class SexsPtiInputs : size_t
    {
      vt,   ///< \f$V_t\f$ Terminal voltage, line-to-line RMS [V]
      vref, ///< \f$V_{\mathrm{ref}}\f$ Voltage reference [p.u.]
      vs,   ///< \f$V_S\f$ Stabilizer signal [p.u.]
      vuel, ///< \f$V_{\mathrm{uel}}\f$ Underexcitation limiter signal [p.u.]
      voel, ///< \f$V_{\mathrm{oel}}\f$ Overexcitation limiter signal [p.u.]
    };

    /// Parameters of a `SexsPti`
    enum class SexsPtiParameters : size_t
    {
      Ta,     ///< \f$T_A\f$ Numerator time constant of lead-lag block [s]
      Tb,     ///< \f$T_B\f$ Denominator time constant of lead-lag block [s]
      Te,     ///< \f$T_E\f$ Exciter field time constant [s]
      K,      ///< \f$K\f$ Voltage regulator gain [p.u.]
      Efdmax, ///< \f$E_{\mathrm{fd}}^{\max}\f$ Maximum excitation output [p.u.]
      Efdmin, ///< \f$E_{\mathrm{fd}}^{\min}\f$ Minimum excitation output [p.u.]
      V,      ///< \f$V\f$ Rated line-to-line RMS voltage [V]
    };

    /// Data of a `SexsPti`
    template <typename RealT, typename IdxT>
    struct SexsPtiData
    {
      using Parameters = SexsPtiParameters;

      RealT Ta{ZERO<RealT>};
      RealT Tb{ZERO<RealT>};
      RealT Te{ZERO<RealT>};
      RealT K{ZERO<RealT>};
      RealT Efdmax{ZERO<RealT>};
      RealT Efdmin{ZERO<RealT>};
      RealT V{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::Ta, &SexsPtiData::Ta},
                                             Parameter{Parameters::Tb, &SexsPtiData::Tb},
                                             Parameter{Parameters::Te, &SexsPtiData::Te},
                                             Parameter{Parameters::K, &SexsPtiData::K},
                                             Parameter{Parameters::Efdmax, &SexsPtiData::Efdmax},
                                             Parameter{Parameters::Efdmin, &SexsPtiData::Efdmin},
                                             Parameter{Parameters::V, &SexsPtiData::V}};
    };
  } // namespace EMT
} // namespace GridKit
