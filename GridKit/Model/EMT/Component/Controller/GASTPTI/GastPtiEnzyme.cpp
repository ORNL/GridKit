/**
 * @file GastPtiEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include <array>
#include <map>

#include <GridKit/Model/EMT/SignalJacobian.hpp>

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /**
       * @brief Sparse Jacobian of the internal equations and attached signals
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::evaluateJacobian()
      {
        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks.
          // The size of the buffer is the sum of maximum capacities of the blocks.
          // Enzyme will compute the appropriate nnz from sparsification.
          auto size         = static_cast<size_t>(size_);
          auto y_ext_size   = y_ext_.size();
          auto buffer_size  = 2 * size * size + size * y_ext_size;
          buffer_size      *= this->externalJacobianExpansion();
          J_rows_buffer_    = new IdxT[buffer_size];
          J_cols_buffer_    = new IdxT[buffer_size];
          J_vals_buffer_    = new RealT[buffer_size];
        }

        nnz_ = 0;

        using GridKit::Enzyme::Sparse::Equation;
        using GridKit::Enzyme::Sparse::SparseJacobian;
        using GridKit::Enzyme::Sparse::Variable;
        using ModelT = GastPti<ScalarT, IdxT>;

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

        SparseJacobian<ModelT, Equation::Internal, Variable::Y>::eval(
            this, n_f, n_y, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);
        SparseJacobian<ModelT, Equation::Internal, Variable::Yp>::eval(
            this, n_f, n_y, ri, vi, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_, alpha_);
        SignalJacobian<ModelT, Equation::Internal, Variable::YExt>::eval(
            this, this, n_f, n_yext, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);

        // Keep structural zeros when a limiter changes its active region. System assembly reuses this COO layout.
        std::map<std::pair<IdxT, IdxT>, RealT>          entries;
        constexpr std::array<std::pair<IdxT, IdxT>, 14> internal_pattern{{{0, 0}, {0, 5}, {1, 0}, {1, 1}, {2, 1}, {2, 2}, {3, 3}, {4, 2}, {4, 4}, {5, 3}, {5, 4}, {5, 5}, {6, 1}, {6, 6}}};
        for (const auto& [row, column] : internal_pattern)
          entries[{ri[row], vi[column]}] = ZERO<RealT>;
        constexpr std::array<std::pair<IdxT, size_t>, 3> external_pattern{{{3, 0}, {3, 1}, {6, 0}}};
        for (const auto& [row, slot] : external_pattern)
        {
          auto* signal = this->externalVariableSignals()[slot];
          if (signal != nullptr)
          {
            typename SignalT::GradientT gradient;
            signal->appendGradient(gradient);
            for (const auto& [column, coefficient] : gradient)
            {
              (void) coefficient;
              entries[{ri[row], column}] = ZERO<RealT>;
            }
          }
        }
        for (IdxT i = 0; i < nnz_; ++i)
          entries[{J_rows_buffer_[i], J_cols_buffer_[i]}] += J_vals_buffer_[i];
        nnz_ = 0;
        for (const auto& [position, value] : entries)
        {
          J_rows_buffer_[nnz_]   = position.first;
          J_cols_buffer_[nnz_]   = position.second;
          J_vals_buffer_[nnz_++] = value;
        }
        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class GastPti<double, long int>;
      template class GastPti<double, size_t>;

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
