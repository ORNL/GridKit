/**
 * @file ParkData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT Park transforms.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `Park`
    enum class ParkInputs : size_t
    {
      abc,   ///< \f$\mathbf{u}\f$ Three-phase quantities, phases a, b and c
      theta, ///< \f$\theta\f$ Angle of the d axis from the axis of phase a [rad]
    };

    /// Parameters of a `Park`
    enum class ParkParameters : size_t
    {
      base, ///< \f$U^\mathrm{base}\f$ Rated peak magnitude, in the units of the input
    };

    /// Data of a `Park`
    template <typename RealT, typename IdxT>
    struct ParkData
    {
      using Parameters = ParkParameters;

      RealT base{ONE<RealT>};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::base, &ParkData::base}};
    };
  } // namespace EMT
} // namespace GridKit
