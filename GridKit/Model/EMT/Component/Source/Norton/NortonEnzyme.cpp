/**
 * @file NortonEnzyme.cpp
 * @brief Enzyme Jacobian of the EMT Norton model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobian.hpp>

#include "NortonImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      const size_t capacity = 9 + static_cast<size_t>(admittance_.jacobianCapacity()) * admittance_.externalJacobianExpansion();
      if (capacity > jacobian_capacity_)
      {
        this->resetJacobianStructure();
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
      using ModelT = Norton<ScalarT, IdxT>;
      if (y_scale != ZERO<RealT>)
      {
        SparseJacobian<ModelT, Equation::Internal, Variable::Y>::eval(
            this, 3, 3, residual_indices_.data(), variable_indices_.data(), y_.getData(), nullptr, nullptr, nullptr, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
      }
      const int status = this->evaluateOperatorJacobians(y_scale, yp_scale);
      if (status != 0)
        return status;
      this->appendOperatorJacobians();
      return this->constructCoo();
    }

    template class Norton<double, long int>;
    template class Norton<double, size_t>;
  } // namespace EMT
} // namespace GridKit
