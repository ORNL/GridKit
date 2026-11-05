/**
 * @file PmuData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT phasor measurement units.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pll/PllData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Pmu`
    enum class PmuInputs : size_t
    {
      v,   ///< \f$\mathbf{v}\f$ Voltage that is measured, phases a, b and c [V]
      i,   ///< \f$\mathbf{i}\f$ Current that is measured, phases a, b and c [A]
      pll, ///< Angle and integrator state of submodel `pll`, which the unit connects
    };

    /// Parameters of a `Pmu`
    enum class PmuParameters : size_t
    {
      pll, ///< Phase-locked loop
    };

    /// Data of a `Pmu`
    template <typename RealT, typename IdxT>
    struct PmuData
    {
      using Parameters = PmuParameters;

      PllData<RealT, IdxT> pll;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::pll, &PmuData::pll}};
    };
  } // namespace EMT
} // namespace GridKit
