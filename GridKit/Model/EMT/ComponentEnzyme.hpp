/**
 * @file ComponentEnzyme.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme partials of fixed-layout EMT models.
 *
 * Enzyme's sparsity analysis needs the local values on the stack and the
 * numbers of internal and external variables as run-time loop bounds.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <new>

#include <GridKit/AutomaticDifferentiation/Enzyme/DfDy.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::tagDifferentiable(const ModelT& model)
    {
      const int status = evaluateDfDy(model, ZERO<RealT>, ONE<RealT>);
      if (dfdyp_.resize(jacobian_entries_.nnz) != 0)
      {
        Log::error() << id_ << ": failed to allocate the derivative entries\n";
        throw std::bad_alloc{};
      }
      std::copy_n(jacobian_entries_.values.getData(), jacobian_entries_.nnz, dfdyp_.getData());
      return status;
    }

    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateJacobian(const ModelT& model)
    {
      return evaluateDfDy(model, ONE<RealT>, alpha_);
    }

    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateDfDy(const ModelT& model, RealT beta, RealT alpha)
    {
      using DfDyT = Enzyme::Sparse::DfDy<ModelT, Enzyme::Sparse::MemberFunctions::Residual>;

      std::array<ScalarT, localSize<ModelT>()> y;
      std::array<ScalarT, localSize<ModelT>()> yp;
      gather(y.data(), yp.data());

      jacobian_entries_.nnz = 0;
      DfDyT::eval(&model,
                  size_,
                  variable_indices_.getSize() - size_,
                  variable_indices_.getData(),
                  y.data(),
                  yp.data(),
                  beta,
                  alpha,
                  jacobian_entries_.rows.getData(),
                  jacobian_entries_.cols.getData(),
                  jacobian_entries_.values.getData(),
                  jacobian_entries_.nnz);
      assert(jacobian_entries_.nnz <= jacobian_entries_.values.getSize());
      return 0;
    }
  } // namespace EMT
} // namespace GridKit
