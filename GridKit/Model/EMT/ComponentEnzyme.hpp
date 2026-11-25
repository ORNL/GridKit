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
    /**
     * @brief Store dF/dy' of `kernels` at the current values.
     *
     * @tparam KernelsT Fixed-layout residual kernels: a model, or a set that a model selects among.
     * @param[in] kernels Kernels to differentiate.
     * @return 0.
     *
     * @throws std::bad_alloc if the partials cannot be allocated.
     *
     * @pre The system has bound this component and resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> KernelsT>
    int Component<scalar_type, index_type>::tagDifferentiable(const KernelsT& kernels)
    {
      const int status = evaluateDfDy(kernels, ZERO<RealT>, ONE<RealT>);
      if (dfdyp_.resize(jacobian_entries_.nnz) != 0)
      {
        Log::error() << id_ << ": failed to allocate the derivative entries\n";
        throw std::bad_alloc{};
      }
      std::copy_n(jacobian_entries_.values.getData(), jacobian_entries_.nnz, dfdyp_.getData());
      return status;
    }

    /**
     * @brief Set the Jacobian entries of `kernels` to dF/dy + alpha dF/dy' at
     *        the current values.
     *
     * @tparam KernelsT Fixed-layout residual kernels: a model, or a set that a model selects among.
     * @param[in] kernels Kernels to differentiate.
     * @return 0.
     *
     * @pre The system has bound this component and resolved its indices.
     */
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> KernelsT>
    int Component<scalar_type, index_type>::evaluateJacobian(const KernelsT& kernels)
    {
      return evaluateDfDy(kernels, ONE<RealT>, clock_->alpha);
    }

    /**
     * @brief Set the Jacobian entries of `model` to beta dF/dy + alpha dF/dy'
     *        at the current values.
     *
     * The local values are copied to the stack: the internal variables, then
     * the external values of the terminals and of the inputs. An output that
     * is not connected is not differentiated and has no entries.
     *
     * @tparam ModelT Fixed-layout model providing the residual kernels.
     * @param[in] model This component as its derived type.
     * @param[in] beta Coefficient of dF/dy.
     * @param[in] alpha Coefficient of dF/dy'.
     * @return 0.
     *
     * @pre The system has bound this component and resolved its indices.
     * @post `jacobian_entries_` holds the entries in system indices.
     */
    template <typename scalar_type, typename index_type>
    template <ComponentKernels<scalar_type, index_type> ModelT>
    int Component<scalar_type, index_type>::evaluateDfDy(const ModelT& model, RealT beta, RealT alpha)
    {
      using DfDyT = Enzyme::Sparse::DfDy<ModelT, Enzyme::Sparse::MemberFunctions::Residual>;

      std::array<ScalarT, valueSize<ModelT>()> y;
      std::array<ScalarT, valueSize<ModelT>()> yp;
      const auto [w, wp] = gather();
      std::copy_n(y_.getData(), size_, y.data());
      std::copy_n(yp_.getData(), size_, yp.data());
      std::copy_n(w, external_indices_.size(), y.data() + size_);
      std::copy_n(wp, external_indices_.size(), yp.data() + size_);

      IdxT*  rows   = jacobian_entries_.rows.getData();
      IdxT*  cols   = jacobian_entries_.cols.getData();
      RealT* values = jacobian_entries_.values.getData();

      jacobian_entries_.nnz = 0;
      DfDyT::eval(&model,
                  size_,
                  external_indices_.size(),
                  residual_indices_.data(),
                  variable_indices_.getData(),
                  y.data(),
                  yp.data(),
                  beta,
                  alpha,
                  rows,
                  cols,
                  values,
                  jacobian_entries_.nnz);
      assert(jacobian_entries_.nnz <= jacobian_entries_.values.getSize());
      return 0;
    }
  } // namespace EMT
} // namespace GridKit
