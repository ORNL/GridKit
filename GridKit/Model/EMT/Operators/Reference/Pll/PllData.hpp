/**
 * @file PllData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT phase-locked loops.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Pll`
    enum class PllInputs : size_t
    {
      v, ///< \f$\mathbf{v}\f$ Voltage that is tracked, phases a, b and c [V]
    };

    /// Parameters of a `Pll`
    enum class PllParameters : size_t
    {
      Kp, ///< \f$K_P\f$ Proportional gain [rad/(V s)]
      Ki, ///< \f$K_I\f$ Integral gain [rad/(V s^2)]
    };

    /// Data of a `Pll`
    template <typename RealT, typename IdxT>
    struct PllData
    {
      using Parameters = PllParameters;

      RealT Kp{ZERO<RealT>};
      RealT Ki{ZERO<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::Kp, &PllData::Kp}, Parameter{Parameters::Ki, &PllData::Ki}};
    };
  } // namespace EMT
} // namespace GridKit
