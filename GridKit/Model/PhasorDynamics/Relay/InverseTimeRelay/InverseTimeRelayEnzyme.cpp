/**
 * @file InverseTimeRelayEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the InverseTimeRelay model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "InverseTimeRelayImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Relay
    {
      /**
       * @brief Assemble the sparse InverseTimeRelay Jacobian with Enzyme.
       *
       * Differentiates the internal residual with respect to the internal
       * variables, their derivatives, and the measured current. A measured
       * current read from a computed output reaches the system Jacobian
       * through the chain rule in SystemModel.
       *
       * @pre allocate() has completed.
       * @pre evaluateResidual() has refreshed the signal buffer at the
       *      current state.
       * @pre Solver alpha and global variable and residual indices are set.
       */
      template <typename scalar_type, typename index_type>
      int InverseTimeRelay<scalar_type, index_type>::evaluateJacobian()
      {
        using ModelT = InverseTimeRelay<ScalarT, IdxT>;
        using Fn     = Enzyme::Sparse::MemberFunctions;

        const auto n_y  = static_cast<size_t>(size_);
        const auto n_ws = static_cast<size_t>(ws_.getSize());

        if (J_rows_buffer_ == nullptr)
        {
          const auto buffer_size = 2 * n_y * n_y + n_y * n_ws;
          J_rows_buffer_         = new IdxT[buffer_size];
          J_cols_buffer_         = new IdxT[buffer_size];
          J_vals_buffer_         = new RealT[buffer_size];
        }

        ScalarT* y  = y_.getData();
        ScalarT* yp = yp_.getData();
        ScalarT* ws = ws_.getData();

        nnz_ = 0;

        Enzyme::Sparse::DfDy<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                           n_y,
                                                                           n_y,
                                                                           residual_indices_.data(),
                                                                           variable_indices_.data(),
                                                                           y,
                                                                           yp,
                                                                           nullptr,
                                                                           ws,
                                                                           J_rows_buffer_,
                                                                           J_cols_buffer_,
                                                                           J_vals_buffer_,
                                                                           nnz_);

        Enzyme::Sparse::DfDyp<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                            n_y,
                                                                            n_y,
                                                                            residual_indices_.data(),
                                                                            variable_indices_.data(),
                                                                            y,
                                                                            yp,
                                                                            nullptr,
                                                                            ws,
                                                                            alpha_,
                                                                            J_rows_buffer_,
                                                                            J_cols_buffer_,
                                                                            J_vals_buffer_,
                                                                            nnz_);

        Enzyme::Sparse::DfDws<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                            n_y,
                                                                            n_ws,
                                                                            residual_indices_.data(),
                                                                            ws_indices_.data(),
                                                                            y,
                                                                            yp,
                                                                            nullptr,
                                                                            ws,
                                                                            J_rows_buffer_,
                                                                            J_cols_buffer_,
                                                                            J_vals_buffer_,
                                                                            nnz_);

        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class InverseTimeRelay<double, long int>;
      template class InverseTimeRelay<double, size_t>;
    } // namespace Relay
  } // namespace PhasorDynamics
} // namespace GridKit
