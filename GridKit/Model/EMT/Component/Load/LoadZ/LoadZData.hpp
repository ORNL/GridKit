/**
 * @file LoadZData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT impedance loads.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Ports of a `LoadZ`, each of size N
    enum class LoadZPorts : size_t
    {
      bus, ///< Bus of the load
    };

    /// Parameters of a `LoadZ`
    enum class LoadZParameters : size_t
    {
      N, ///< \f$N\f$ Number of phases
      Z, ///< \f$\mathbf{z}\f$ Impedance, N by N
    };

    /// Data of a `LoadZ`
    template <typename RealT, typename IdxT>
    struct LoadZData
    {
      using Parameters = LoadZParameters;

      IdxT                       N{0};
      VectorFitData<RealT, IdxT> Z;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &LoadZData::N}, Parameter{Parameters::Z, &LoadZData::Z}};
    };
  } // namespace EMT
} // namespace GridKit
