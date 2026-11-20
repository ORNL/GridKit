/**
 * @file VectorFitData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT VectorFit operators.
 */

#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Inputs of a `VectorFit`, connected by its parent
    enum class VectorFitInputs : size_t
    {
      u, ///< \f$\mathbf{u}\f$ Input, of size K
    };

    /// Outputs of a `VectorFit`, connected by its parent
    enum class VectorFitOutputs : size_t
    {
      y, ///< \f$\mathbf{y}\f$ Output, of size N
    };

    /// Parameters of a `VectorFit` with output dimension N and input dimension K
    enum class VectorFitParameters : size_t
    {
      D,        ///< \f$\mathbf{D}\f$ Constant coefficient, N by K [[y]/[u]]
      E,        ///< \f$\mathbf{E}\f$ Linear coefficient, N by K [s [y]/[u]]
      poles,    ///< \f$\mathbf{p}\f$ Poles [1/s]
      residues, ///< \f$\mathbf{R}\f$ Residues, N by K by Q [[y]/(s [u])]
    };

    /// Data of a `VectorFit`
    template <typename RealT, typename IdxT>
    struct VectorFitData
    {
      using Parameters = VectorFitParameters;

      Matrix<RealT>                    D;
      Matrix<RealT>                    E;
      std::vector<Pair<RealT>>         poles;
      Matrix<std::vector<Pair<RealT>>> residues;

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::D, &VectorFitData::D},
                                             Parameter{Parameters::E, &VectorFitData::E},
                                             Parameter{Parameters::poles, &VectorFitData::poles},
                                             Parameter{Parameters::residues, &VectorFitData::residues}};
    };
  } // namespace EMT
} // namespace GridKit
