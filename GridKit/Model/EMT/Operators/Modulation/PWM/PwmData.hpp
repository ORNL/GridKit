/**
 * @file PwmData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT pulse-width modulation.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Parameters of a `Pwm`
    enum class PwmParameters : size_t
    {
      N,         ///< \f$N\f$ Number of channels
      fc,        ///< \f$f_\mathrm{c}\f$ Carrier frequency [Hz]
      alignment, ///< \f$\alpha\f$ Pulse alignment
      injection, ///< \f$\mathrm{injection}\f$ Zero-sequence injection
      k0,        ///< \f$k_0\f$ Zero-state partition
      mu,        ///< \f$\mu\f$ Switching-edge smoothing scale [1/s]
    };

    /// Data of a `Pwm`
    template <typename RealT, typename IdxT>
    struct PwmData
    {
      using Parameters = PwmParameters;

      IdxT  N{0};
      RealT fc{ZERO<RealT>};
      RealT alignment{HALF<RealT>};
      bool  injection{true};
      RealT k0{HALF<RealT>};
      RealT mu{240.0};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &PwmData::N},
                                             Parameter{Parameters::fc, &PwmData::fc},
                                             Parameter{Parameters::alignment, &PwmData::alignment},
                                             Parameter{Parameters::injection, &PwmData::injection},
                                             Parameter{Parameters::k0, &PwmData::k0},
                                             Parameter{Parameters::mu, &PwmData::mu}};
    };
  } // namespace EMT
} // namespace GridKit
