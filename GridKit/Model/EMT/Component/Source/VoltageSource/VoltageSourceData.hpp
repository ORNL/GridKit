/**
 * @file VoltageSourceData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT voltage sources.
 */

#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

#include <GridKit/Model/EMT/ComponentData.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Ports of a `VoltageSource`, each of size N
    enum class VoltageSourcePorts : size_t
    {
      bus, ///< Bus of the source
    };

    /// Parameters of a `VoltageSource`
    enum class VoltageSourceParameters : size_t
    {
      N,     ///< \f$N\f$ Number of phases
      E,     ///< \f$\mathbf{E}\f$ RMS source voltages [V]
      phi,   ///< \f$\boldsymbol{\phi}\f$ Source phase offsets [rad]
      omega, ///< \f$\omega\f$ Source angular frequency [rad/s]
      Y,     ///< \f$\mathbf{y}\f$ Terminal admittance, N by N
    };

    /// Data of a `VoltageSource`
    template <typename RealT, typename IdxT>
    struct VoltageSourceData
    {
      using Parameters = VoltageSourceParameters;

      IdxT                       N{0};
      std::vector<RealT>         E;
      std::vector<RealT>         phi;
      RealT                      omega{ZERO<RealT>};
      VectorFitData<RealT, IdxT> Y;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &VoltageSourceData::N},
                                             Parameter{Parameters::E, &VoltageSourceData::E},
                                             Parameter{Parameters::phi, &VoltageSourceData::phi},
                                             Parameter{Parameters::omega, &VoltageSourceData::omega},
                                             Parameter{Parameters::Y, &VoltageSourceData::Y}};
    };
  } // namespace EMT
} // namespace GridKit
