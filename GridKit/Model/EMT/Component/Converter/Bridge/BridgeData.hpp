/**
 * @file BridgeData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT converter bridges.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/PWM/PwmData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Terminals of a `Bridge`
    enum class BridgeTerminals : size_t
    {
      bus,  ///< Bus of the grid-side filter, phases a, b and c
      vpos, ///< Positive DC terminal
      vneg, ///< Negative DC terminal
    };

    /// Inputs of a `Bridge`
    enum class BridgeInputs : size_t
    {
      theta, ///< \f$\theta\f$ Angle of the dq coordinates [rad]
      mdq,   ///< \f$\mathbf{m}_\mathrm{dq}\f$ Modulation signals, d then q
    };

    /// Parameters of a `Bridge`
    enum class BridgeParameters : size_t
    {
      N,   ///< \f$N\f$ Number of legs
      Rs,  ///< \f$R_\mathrm{s}\f$ Converter-side resistance [Ohm]
      Ls,  ///< \f$L_\mathrm{s}\f$ Converter-side inductance [H]
      C,   ///< \f$C\f$ Filter capacitance [F]
      Rg,  ///< \f$R_g\f$ Grid-side resistance [Ohm]
      Lg,  ///< \f$L_g\f$ Grid-side inductance [H]
      pwm, ///< \f$\mathbf{p}\f$ Pulse-width modulation
    };

    /// Data of a `Bridge`
    template <typename RealT, typename IdxT>
    struct BridgeData
    {
      using Parameters = BridgeParameters;

      IdxT                 N{0};
      RealT                Rs{ZERO<RealT>};
      RealT                Ls{ZERO<RealT>};
      RealT                C{ZERO<RealT>};
      RealT                Rg{ZERO<RealT>};
      RealT                Lg{ZERO<RealT>};
      PwmData<RealT, IdxT> pwm;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &BridgeData::N},
                                             Parameter{Parameters::Rs, &BridgeData::Rs},
                                             Parameter{Parameters::Ls, &BridgeData::Ls},
                                             Parameter{Parameters::C, &BridgeData::C},
                                             Parameter{Parameters::Rg, &BridgeData::Rg},
                                             Parameter{Parameters::Lg, &BridgeData::Lg},
                                             Parameter{Parameters::pwm, &BridgeData::pwm}};
    };
  } // namespace EMT
} // namespace GridKit
