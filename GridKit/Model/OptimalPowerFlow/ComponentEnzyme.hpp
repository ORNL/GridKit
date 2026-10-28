/**
 * @file ComponentEnzyme.hpp
 * @brief Definition of the `Component` Enzyme derivatives.
 *
 * Local variables live on the stack, where Enzyme's sparsity analysis works,
 * and their count is a run-time argument, which Enzyme needs for the loop.
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
        gather(x.data());

        jacobian_entries_.nnz = 0;
        Enzyme::Sparse::DgDx<ModelT>::eval(&model,
                                           local_x_.getSize(),
                                           constraint_indices_.getData(),
                                           variable_indices_.getData(),
                                           x.data(),
                                           jacobian_entries_.rows.getData(),
                                           jacobian_entries_.cols.getData(),
                                           jacobian_entries_.values.getData(),
                                           jacobian_entries_.nnz);
        assert(jacobian_entries_.nnz <= jacobian_entries_.values.getSize());
        return 0;
      }
    }

    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateHessian(const ModelT& model, RealT sigma)
    {
      std::array<ScalarT, localSize<ModelT>().variables> x;
      std::array<RealT, localSize<ModelT>().constraints> lambda;
      gather(x.data());
      gatherMultipliers(lambda.data());

      hessian_entries_.nnz = 0;
      Enzyme::Sparse::D2LDx2<ModelT>::eval(&model,
                                           local_x_.getSize(),
                                           variable_indices_.getData(),
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

    /// Keep the entries on or below the diagonal, which the others mirror
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
