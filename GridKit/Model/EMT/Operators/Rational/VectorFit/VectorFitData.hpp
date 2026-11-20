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

    /// `data` with `D`, `E` and the residues multiplied by `factor`
    template <typename RealT, typename IdxT>
    VectorFitData<RealT, IdxT> scaled(VectorFitData<RealT, IdxT> data, RealT factor)
    {
      for (Matrix<RealT>* matrix : {&data.D, &data.E})
      {
        for (auto& row : *matrix)
        {
          for (RealT& value : row)
          {
            value *= factor;
          }
        }
      }
      for (auto& row : data.residues)
      {
        for (auto& entry : row)
        {
          for (Pair<RealT>& residue : entry)
          {
            residue[0] *= factor;
            residue[1] *= factor;
          }
        }
      }
      return data;
    }

    /// True if the case gives no `D`, `E` or residues
    template <typename RealT, typename IdxT>
    bool omitted(const VectorFitData<RealT, IdxT>& data)
    {
      return data.D.empty() && data.E.empty() && data.residues.empty();
    }
  } // namespace EMT
} // namespace GridKit
