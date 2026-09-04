/**
 * @file Tgov1Enzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobian.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /**
       * @brief Jacobian evaluation experimental
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Tgov1..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        gatherExternalVariables();

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks.
          // The size of the buffer is the sum of maximum capacities of the blocks.
          // Enyme will compute the appropriate nnz from sparsification.
          auto size        = static_cast<size_t>(size_);
          auto y_ext_size  = y_ext_.size();
          auto buffer_size = 2 * size * size + size * y_ext_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        nnz_ = 0;

        using GridKit::Enzyme::Sparse::Equation;
        using GridKit::Enzyme::Sparse::SparseJacobian;
        using GridKit::Enzyme::Sparse::Variable;
        using ModelT = Tgov1<ScalarT, IdxT>;

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
        SparseJacobian<ModelT, Equation::Internal, Variable::YExt>::eval(
            this, n_f, n_yext, ri, vie, y, yp, ye, ype, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz_);

        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
