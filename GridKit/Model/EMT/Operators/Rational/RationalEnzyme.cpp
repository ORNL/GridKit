/**
 * @file RationalEnzyme.cpp
 * @brief Enzyme Jacobian of the EMT Rational model.
 */

#include <GridKit/Model/EMT/SignalJacobian.hpp>

#include "RationalImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      if (!coupling_allocated_ || errors_)
        return 1;
      this->gatherExternalVariables();

      const size_t n_y      = static_cast<size_t>(size_);
      const size_t capacity = 2 * (n_y + rows_) * (n_y + cols_) * this->externalJacobianExpansion();
      if (this->hasComputedInputs() || capacity > capacity_)
        this->resetJacobianStructure();
      if (capacity > capacity_)
      {
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_ = new IdxT[capacity];
        J_cols_buffer_ = new IdxT[capacity];
        J_vals_buffer_ = new RealT[capacity];
        capacity_      = capacity;
      }
      nnz_ = 0;

      using GridKit::Enzyme::Sparse::Equation;
      using GridKit::Enzyme::Sparse::SparseJacobian;
      using GridKit::Enzyme::Sparse::Variable;
      using ModelT = Rational<ScalarT, IdxT>;

      const auto* ri  = residual_indices_.data();
      const auto* rie = residual_indices_ext_.data();
      const auto* vi  = variable_indices_.data();
      const auto* vie = variable_indices_ext_.data();
      const auto* y   = n_y == 0 ? nullptr : y_.getData();
      const auto* yp  = n_y == 0 ? nullptr : yp_.getData();
      const auto* ye  = y_ext_.data();
      const auto* ype = yp_ext_.data();

      if (y_scale != ZERO<RealT>)
      {
        SparseJacobian<ModelT, Equation::Internal, Variable::Y>::eval(
            this, n_y, n_y, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
        SignalJacobian<ModelT, Equation::Internal, Variable::YExt>::eval(
            this, this, n_y, cols_, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
        SparseJacobian<ModelT, Equation::External, Variable::Y>::eval(
            this, rows_, n_y, rie, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
        SignalJacobian<ModelT, Equation::External, Variable::YExt>::eval(
            this, this, rows_, cols_, rie, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, y_scale);
      }

      if (yp_scale != ZERO<RealT>)
      {
        SparseJacobian<ModelT, Equation::Internal, Variable::Yp>::eval(
            this, n_y, n_y, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, yp_scale);
        SignalJacobian<ModelT, Equation::Internal, Variable::YpExt>::eval(
            this, this, n_y, cols_, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, yp_scale);
        SparseJacobian<ModelT, Equation::External, Variable::Yp>::eval(
            this, rows_, n_y, rie, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, yp_scale);
        SignalJacobian<ModelT, Equation::External, Variable::YpExt>::eval(
            this, this, rows_, cols_, rie, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, yp_scale);
      }
      IdxT retained = 0;
      for (IdxT j = 0; j < nnz_; ++j)
      {
        if (J_rows_buffer_[j] == INVALID_INDEX<IdxT>)
          continue;
        J_rows_buffer_[retained]   = J_rows_buffer_[j];
        J_cols_buffer_[retained]   = J_cols_buffer_[j];
        J_vals_buffer_[retained++] = J_vals_buffer_[j];
      }
      nnz_ = retained;
      return this->constructCoo();
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::appendOutputGradient(
        IdxT n, typename SignalT::GradientT& gradient, RealT scale) const
    {
      for (size_t k = 0; k < cols_; ++k)
        if (E_[static_cast<size_t>(n)][k] != ZERO<RealT>)
          throw std::logic_error("A derivative-dependent rational output cannot be a computed signal");
      const size_t n_y = static_cast<size_t>(size_);
      if (output_partials_.empty())
        output_partials_.resize(rows_);
      auto& partials = output_partials_[static_cast<size_t>(n)];
      if (partials.empty())
      {
        auto* model = const_cast<Rational*>(this);
        model->gatherExternalVariables();
        partials.resize(n_y + cols_);
        output_direction_.assign(std::max(n_y, cols_), ScalarT{0});
        const auto* y        = n_y == 0 ? nullptr : y_.getData();
        const auto* ye       = y_ext_.data();
        const auto* ype      = yp_ext_.data();
        auto        evaluate = +[](const Rational* model, IdxT row, const ScalarT* y, const ScalarT* ye, const ScalarT* ype)
        {
          return model->evaluateOutput(row, y, ye, ype);
        };
        for (size_t k = 0; k < n_y; ++k)
        {
          output_direction_[k] = ScalarT{1};
          partials[k]          = Enzyme::Sparse::__enzyme_fwddiff<ScalarT>(
              (void*) evaluate, enzyme_const, this, enzyme_const, n, enzyme_dup, y, output_direction_.data(), enzyme_const, ye, enzyme_const, ype);
          output_direction_[k] = ScalarT{0};
        }
        for (size_t k = 0; k < cols_; ++k)
        {
          output_direction_[k] = ScalarT{1};
          partials[n_y + k]    = Enzyme::Sparse::__enzyme_fwddiff<ScalarT>(
              (void*) evaluate, enzyme_const, this, enzyme_const, n, enzyme_const, y, enzyme_dup, ye, output_direction_.data(), enzyme_const, ype);
          output_direction_[k] = ScalarT{0};
        }
      }
      for (size_t k = 0; k < n_y; ++k)
        if (partials[k] != ZERO<RealT>)
          gradient.emplace_back(this->getVariableIndex(static_cast<IdxT>(k)), scale * partials[k]);
      for (size_t k = 0; k < cols_; ++k)
        if (partials[n_y + k] != ZERO<RealT>)
          input_[k]->appendGradient(gradient, scale * partials[n_y + k]);
    }

    template class Rational<double, long int>;
    template class Rational<double, size_t>;
  } // namespace EMT
} // namespace GridKit
