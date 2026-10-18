/**
 * @file ComponentEnzyme.hpp
 * @brief Definition of the `Component` Enzyme derivative evaluations.
 *
 * The derivatives gather the local variables onto the stack, where Enzyme's
 * sparsity analysis is more accurate than over the heap buffer, and pass
 * their count at run time, which Enzyme needs to find the derivative loop.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>

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
    int Component<scalar_type, index_type>::evaluateEqualityJacobian(const ModelT& model)
    {
      if constexpr (localSize<ModelT>().equalities == 0)
      {
        equality_entries_.nnz = 0;
        return 0;
      }
      else
      {
        using DriverT = Enzyme::Sparse::DgDx<ModelT, Enzyme::Sparse::MemberFunctions::Equalities>;

        std::array<IdxT, localSize<ModelT>().equalities> rows;
        equalityIndices(rows.data());
        return evaluateJacobian<DriverT>(model, rows.data(), equality_entries_);
      }
    }

    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateInequalityJacobian(const ModelT& model)
    {
      if constexpr (localSize<ModelT>().inequalities == 0)
      {
        inequality_entries_.nnz = 0;
        return 0;
      }
      else
      {
        using DriverT = Enzyme::Sparse::DgDx<ModelT, Enzyme::Sparse::MemberFunctions::Inequalities>;

        std::array<IdxT, localSize<ModelT>().inequalities> rows;
        std::iota(rows.begin(), rows.end(), offset_.inequalities);
        return evaluateJacobian<DriverT>(model, rows.data(), inequality_entries_);
      }
    }

    /**
     * @brief Lower triangle of the Hessian of the component Lagrangian
     *
     * @param[in] lambda_c - Equality multipliers in global rows
     * @param[in] lambda_d - Inequality multipliers in global rows, contiguous
     * for the internal block
     */
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateHessian(const ModelT& model,
                                                            RealT         sigma,
                                                            const RealT*  lambda_c,
                                                            const RealT*  lambda_d)
    {
      std::array<ScalarT, localSize<ModelT>().variables> x;
      std::array<IdxT, localSize<ModelT>().variables>    cols;
      std::array<IdxT, localSize<ModelT>().equalities>   rows;
      std::array<RealT, localSize<ModelT>().equalities>  multipliers;
      gather(x.data());
      variableIndices(cols.data());
      equalityIndices(rows.data());
      std::ranges::transform(rows, multipliers.begin(), [lambda_c](IdxT row)
                             { return lambda_c[row]; });

      hessian_entries_.nnz = 0;
      Enzyme::Sparse::D2LDx2<ModelT>::eval(&model,
                                           local_x_.getSize(),
                                           cols.data(),
                                           x.data(),
                                           sigma,
                                           multipliers.data(),
                                           lambda_d + offset_.inequalities,
                                           hessian_entries_.rows.getData(),
                                           hessian_entries_.cols.getData(),
                                           hessian_entries_.values.getData(),
                                           hessian_entries_.nnz);
      assert(hessian_entries_.nnz <= hessian_entries_.values.getSize());

      keepLowerTriangle(hessian_entries_);
      return 0;
    }

    /**
     * @brief Jacobian entries of the kernel that `DriverT` differentiates, in
     * the global `rows`
     */
    template <typename scalar_type, typename index_type>
    template <typename DriverT, typename ModelT>
    int Component<scalar_type, index_type>::evaluateJacobian(const ModelT& model, const IdxT* rows, CooBuffer& entries)
    {
      std::array<ScalarT, localSize<ModelT>().variables> x;
      std::array<IdxT, localSize<ModelT>().variables>    cols;
      gather(x.data());
      variableIndices(cols.data());

      entries.nnz = 0;
      DriverT::eval(&model,
                    local_x_.getSize(),
                    rows,
                    cols.data(),
                    x.data(),
                    entries.rows.getData(),
                    entries.cols.getData(),
                    entries.values.getData(),
                    entries.nnz);
      assert(entries.nnz <= entries.values.getSize());
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
