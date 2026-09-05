/**
 * @file VectorFitEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include <GridKit/Model/EMT/SignalJacobian.hpp>

#include "VectorFitImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * Every differentiated callee is a fixed-size straight-line section
     * kernel; the runtime pole count lives only in this orchestration loop,
     * which slices the shared index arrays per section. The operator is
     * linear time-invariant, so Enzyme runs only at structure discovery: the
     * alpha-scaled blocks are evaluated last at unit scaling and their range
     * is recorded, and every later call refreshes that range with the
     * current alpha and touches nothing else.
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for VectorFit..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      this->gatherExternalVariables();

      if (jacobian_cached_ && !this->hasComputedInputs())
      {
        for (IdxT i = alpha_range_begin_; i < alpha_range_end_; ++i)
        {
          J_vals_buffer_[i] = alpha_ * unit_vals_[static_cast<size_t>(i - alpha_range_begin_)];
        }
        return 0;
      }

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enyme will compute the appropriate nnz from sparsification.
        auto buffer_size  = static_cast<size_t>(jacobianCapacity());
        buffer_size      *= this->externalJacobianExpansion();
        J_rows_buffer_    = new IdxT[buffer_size];
        J_cols_buffer_    = new IdxT[buffer_size];
        J_vals_buffer_    = new RealT[buffer_size];
      }

      if (this->hasComputedInputs())
      {
        this->resetJacobianStructure();
      }
      nnz_ = 0;

      using GridKit::Enzyme::Sparse::Equation;
      using GridKit::Enzyme::Sparse::SparseJacobian;
      using GridKit::Enzyme::Sparse::Variable;
      using RealSectionT        = VectorFitRealSection<ScalarT, IdxT>;
      using ComplexSectionT     = VectorFitComplexSection<ScalarT, IdxT>;
      using FeedthroughSectionT = VectorFitFeedthroughSection<ScalarT, IdxT>;

      const auto* ri_base = (this->getResidualIndices()).data();
      const auto* vi_base = (this->getVariableIndices()).data();
      const auto* rie     = residual_indices_ext_.data();
      const auto* vie     = variable_indices_ext_.data();
      const auto* y_base  = y_.getData();
      const auto* yp_base = yp_.getData();
      const auto* ye      = y_ext_.data();
      const auto* ype     = yp_ext_.data();

      // Constant blocks first: state coupling, input coupling, and output rows
      IdxT offset = 0;
      for (auto& section : real_sections_)
      {
        const auto* ri = ri_base + offset;
        const auto* vi = vi_base + offset;
        const auto* y  = y_base + offset;
        const auto* yp = yp_base + offset;

        SparseJacobian<RealSectionT, Equation::Internal, Variable::Y>::eval(
            &section, 3, 3, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        SignalJacobian<RealSectionT, Equation::Internal, Variable::YExt>::eval(
            this, &section, 3, 3, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        SparseJacobian<RealSectionT, Equation::External, Variable::Y>::eval(
            &section, 3, 3, rie, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        offset += 3;
      }

      for (auto& section : complex_sections_)
      {
        const auto* ri = ri_base + offset;
        const auto* vi = vi_base + offset;
        const auto* y  = y_base + offset;
        const auto* yp = yp_base + offset;

        SparseJacobian<ComplexSectionT, Equation::Internal, Variable::Y>::eval(
            &section, 6, 6, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        SignalJacobian<ComplexSectionT, Equation::Internal, Variable::YExt>::eval(
            this, &section, 6, 3, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        SparseJacobian<ComplexSectionT, Equation::External, Variable::Y>::eval(
            &section, 3, 6, rie, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        offset += 6;
      }

      SignalJacobian<FeedthroughSectionT, Equation::External, Variable::YExt>::eval(
          this, &feedthrough_, 3, 3, rie, vie, y_base, yp_base, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);

      // Alpha-scaled blocks last, contiguous, at unit scaling
      alpha_range_begin_ = nnz_;

      offset = 0;
      for (auto& section : real_sections_)
      {
        SparseJacobian<RealSectionT, Equation::Internal, Variable::Yp>::eval(
            &section, 3, 3, ri_base + offset, vi_base + offset, y_base + offset, yp_base + offset, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, ONE<RealT>);
        offset += 3;
      }
      for (auto& section : complex_sections_)
      {
        SparseJacobian<ComplexSectionT, Equation::Internal, Variable::Yp>::eval(
            &section, 6, 6, ri_base + offset, vi_base + offset, y_base + offset, yp_base + offset, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, ONE<RealT>);
        offset += 6;
      }
      if (hasFeedthroughDerivative())
      {
        SignalJacobian<FeedthroughSectionT, Equation::External, Variable::YpExt>::eval(
            this, &feedthrough_, 3, 3, rie, vie, y_base, yp_base, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, ONE<RealT>);
      }

      alpha_range_end_ = nnz_;

      unit_vals_.assign(J_vals_buffer_ + alpha_range_begin_, J_vals_buffer_ + alpha_range_end_);
      for (IdxT i = alpha_range_begin_; i < alpha_range_end_; ++i)
      {
        J_vals_buffer_[i] = alpha_ * unit_vals_[static_cast<size_t>(i - alpha_range_begin_)];
      }
      jacobian_cached_ = true;

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class VectorFit<double, long int>;
    template class VectorFit<double, size_t>;

  } // namespace EMT
} // namespace GridKit
