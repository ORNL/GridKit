/**
 * @file ComponentEnzyme.hpp
 * @brief Definition of the `Component` Enzyme derivative evaluations.
 *
 * The derivatives gather the local variables onto the stack, where Enzyme's
 * sparsity analysis is more accurate than over the heap buffer, and pass
 * their count at run time, which Enzyme needs to find the derivative loop.
 */

#pragma once

#include <array>
#include <cassert>

#include <GridKit/AutomaticDifferentiation/Enzyme/D2LDx2.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DfDx.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/DgDx.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateGradient(const ModelT& model)
    {
      std::array<ScalarT, localSize<ModelT>().variables> x;
      std::array<RealT, localSize<ModelT>().variables>   gradient;
      gather(x.data());
      Enzyme::Sparse::DfDx<ModelT>::eval(&model, x.data(), gradient.data());

      scatterGradient(gradient.data());
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateJacobian(const ModelT& model)
    {
      if constexpr (localSize<ModelT>().constraints == 0)
      {
        jacobian_entries_.nnz = 0;
        return 0;
      }
      else
      {
        std::array<ScalarT, localSize<ModelT>().variables> x;
        std::array<IdxT, localSize<ModelT>().variables>    cols;
        std::array<IdxT, localSize<ModelT>().constraints>  rows;
        gather(x.data());
        variableIndices(cols.data());
        constraintIndices(rows.data());

        jacobian_entries_.nnz = 0;
        Enzyme::Sparse::DgDx<ModelT>::eval(&model,
                                           local_x_.getSize(),
                                           rows.data(),
                                           cols.data(),
                                           x.data(),
                                           jacobian_entries_.rows.getData(),
                                           jacobian_entries_.cols.getData(),
                                           jacobian_entries_.values.getData(),
                                           jacobian_entries_.nnz);
        assert(jacobian_entries_.nnz <= jacobian_entries_.values.getSize());
        return 0;
      }
    }

    /**
     * @brief Lower triangle of the Hessian of the component Lagrangian at the
     * multipliers of the constraints it contributes to
     */
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateHessian(const ModelT& model, RealT sigma)
    {
      std::array<ScalarT, localSize<ModelT>().variables> x;
      std::array<IdxT, localSize<ModelT>().variables>    cols;
      std::array<RealT, localSize<ModelT>().constraints> lambda;
      gather(x.data());
      variableIndices(cols.data());
      gatherMultipliers(lambda.data());

      hessian_entries_.nnz = 0;
      Enzyme::Sparse::D2LDx2<ModelT>::eval(&model,
                                           local_x_.getSize(),
                                           cols.data(),
                                           x.data(),
                                           sigma,
                                           lambda.data(),
                                           hessian_entries_.rows.getData(),
                                           hessian_entries_.cols.getData(),
                                           hessian_entries_.values.getData(),
                                           hessian_entries_.nnz);
      assert(hessian_entries_.nnz <= hessian_entries_.values.getSize());

      keepLowerTriangle(hessian_entries_);
      return 0;
    }

    /**
     * @brief Keep the entries on or below the diagonal in global indices
     *
     * The symmetric entries above it carry the same values.
     */
    template <typename scalar_type, typename index_type>
    void Component<scalar_type, index_type>::keepLowerTriangle(CooBuffer& entries)
    {
      IdxT*  rows   = entries.rows.getData();
      IdxT*  cols   = entries.cols.getData();
      RealT* values = entries.values.getData();
      IdxT   kept   = 0;
      for (IdxT k = 0; k < entries.nnz; ++k)
      {
        if (rows[k] >= cols[k])
        {
          rows[kept]   = rows[k];
          cols[kept]   = cols[k];
          values[kept] = values[k];
          ++kept;
        }
      }
      entries.nnz = kept;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
