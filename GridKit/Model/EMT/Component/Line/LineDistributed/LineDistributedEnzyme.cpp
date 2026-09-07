/**
 * @file LineDistributedEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include <GridKit/Model/EMT/SignalJacobian.hpp>

#include "LineDistributedImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int LineDistributed<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      this->gatherExternalVariables();

      const int status = this->evaluateOperatorJacobians(y_scale, yp_scale);
      if (status != 0)
        return status;
      const auto size     = static_cast<size_t>(this->equationSize());
      size_t     capacity = 2 * size * (size + y_ext_.size()) * this->externalJacobianExpansion();
      for (auto* op : this->operators_)
        capacity += static_cast<size_t>(op->getCooJacobian()->getNnz());
      if (capacity > jacobian_capacity_)
      {
        delete this->coo_jac_;
        this->coo_jac_ = nullptr;
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_     = new IdxT[capacity];
        J_cols_buffer_     = new IdxT[capacity];
        J_vals_buffer_     = new RealT[capacity];
        jacobian_capacity_ = capacity;
      }

      nnz_ = 0;

      using GridKit::Enzyme::Sparse::Equation;
      using GridKit::Enzyme::Sparse::SparseJacobian;
      using GridKit::Enzyme::Sparse::Variable;
      using ModelT = LineDistributed<ScalarT, IdxT>;

      const auto  n_f    = static_cast<size_t>(f_.getSize());
      const auto  n_y    = static_cast<size_t>(y_.getSize());
      const auto  n_yext = y_ext_.size();
      const auto* ri     = (this->getResidualIndices()).data();
      const auto* vi     = (this->getVariableIndices()).data();
      const auto* vie    = variable_indices_ext_.data();
      const auto* y      = y_.getData();
      const auto* yp     = yp_.getData();
      const auto* ye     = y_ext_.data();
      const auto* ype    = yp_ext_.data();

      if (y_scale != ZERO<RealT>)
      {
        SparseJacobian<ModelT, Equation::Internal, Variable::Y>::eval(
            this, n_f, n_y, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
      }
      if (y_scale != ZERO<RealT>)
      {
        SignalJacobian<ModelT, Equation::Internal, Variable::YExt>::eval(
            this, this, n_f, n_yext, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
      }

      this->appendOperatorJacobians();

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class LineDistributed<double, long int>;
    template class LineDistributed<double, size_t>;

  } // namespace EMT
} // namespace GridKit
