/**
 * @file FilterEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include <GridKit/Model/EMT/SignalJacobian.hpp>

#include "FilterImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Sparse Jacobian of the LCL residual
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Filter<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      this->gatherExternalVariables();

      if (this->J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enzyme will compute the appropriate nnz from sparsification.
        auto size             = static_cast<size_t>(this->size_);
        auto y_ext_size       = this->y_ext_.size();
        auto buffer_size      = 2 * size * (size + y_ext_size);
        buffer_size          *= this->externalJacobianExpansion();
        this->J_rows_buffer_  = new IdxT[buffer_size];
        this->J_cols_buffer_  = new IdxT[buffer_size];
        this->J_vals_buffer_  = new RealT[buffer_size];
      }

      this->nnz_ = 0;

      using GridKit::Enzyme::Sparse::Equation;
      using GridKit::Enzyme::Sparse::SparseJacobian;
      using GridKit::Enzyme::Sparse::Variable;
      using ModelT = Filter<ScalarT, IdxT>;

      const auto  n_f    = static_cast<size_t>(this->f_.getSize());
      const auto  n_y    = static_cast<size_t>(this->y_.getSize());
      const auto  n_yext = this->y_ext_.size();
      const auto* ri     = (this->getResidualIndices()).data();
      const auto* vi     = (this->getVariableIndices()).data();
      const auto* vie    = this->variable_indices_ext_.data();
      const auto* y      = this->y_.getData();
      const auto* yp     = this->yp_.getData();
      const auto* ye     = this->y_ext_.data();
      const auto* ype    = this->yp_ext_.data();

      if (y_scale != ZERO<RealT>)
      {
        SparseJacobian<ModelT, Equation::Internal, Variable::Y>::eval(
            this, n_f, n_y, ri, vi, y, yp, ye, ype, this->J_rows_buffer_, this->J_cols_buffer_, this->J_vals_buffer_, this->nnz_, y_scale);
      }
      if (yp_scale != ZERO<RealT>)
      {
        SparseJacobian<ModelT, Equation::Internal, Variable::Yp>::eval(
            this, n_f, n_y, ri, vi, y, yp, ye, ype, this->J_rows_buffer_, this->J_cols_buffer_, this->J_vals_buffer_, this->nnz_, yp_scale);
      }
      if (y_scale != ZERO<RealT>)
      {
        SignalJacobian<ModelT, Equation::Internal, Variable::YExt>::eval(
            this, this, n_f, n_yext, ri, vie, y, yp, ye, ype, this->J_rows_buffer_, this->J_cols_buffer_, this->J_vals_buffer_, this->nnz_, y_scale);
      }

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class Filter<double, long int>;
    template class Filter<double, size_t>;

  } // namespace EMT
} // namespace GridKit
