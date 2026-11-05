/**
 * @file LineLumpedData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT lumped lines.
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
    /// Ports of a `LineLumped`, each of size N
    enum class LineLumpedPorts : size_t
    {
      bus1, ///< Terminal 1
      bus2, ///< Terminal 2
    };

    /// Parameters of a `LineLumped`
    enum class LineLumpedParameters : size_t
    {
      N,          ///< \f$N\f$ Number of phases
      K,          ///< \f$K\f$ Number of conductors
      conductors, ///< \f$\mathbf{c}\f$ Phase of each conductor, from 1
      dx,         ///< \f$\Delta x\f$ Segment length [m]
      Zp,         ///< \f$\mathbf{z}\f$ Series impedance per unit length, K by K
      Yp,         ///< \f$\mathbf{y}\f$ Shunt admittance per unit length, K by K
    };

    /// Data of a `LineLumped`
    template <typename RealT, typename IdxT>
    struct LineLumpedData
    {
      using Parameters = LineLumpedParameters;

      IdxT                       N{0};
      IdxT                       K{0};
      std::vector<IdxT>          conductors;
      RealT                      dx{ZERO<RealT>};
      VectorFitData<RealT, IdxT> Zp;
      VectorFitData<RealT, IdxT> Yp;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &LineLumpedData::N},
                                             Parameter{Parameters::K, &LineLumpedData::K},
                                             Parameter{Parameters::conductors, &LineLumpedData::conductors},
                                             Parameter{Parameters::dx, &LineLumpedData::dx},
                                             Parameter{Parameters::Zp, &LineLumpedData::Zp},
                                             Parameter{Parameters::Yp, &LineLumpedData::Yp}};
    };
  } // namespace EMT
} // namespace GridKit
