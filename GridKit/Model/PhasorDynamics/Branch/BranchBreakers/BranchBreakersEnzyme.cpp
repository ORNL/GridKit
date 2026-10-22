/**
 * @file BranchBreakersEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the BranchBreakers model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BranchBreakersImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Assemble the sparse BranchBreakers Jacobian with Enzyme.
     *
     * Differentiates the breaker latches, the bus-current blocks with respect
     * to the terminal voltages and the latches, and the terminal-current
     * outputs. SystemModel applies the chain rule to the output rows.
     *
     * @pre allocate() has completed.
     * @pre evaluateResidual() has refreshed the interface buffers at the
     *      current state.
     * @pre Solver alpha and global variable, residual, and output indices
     *      are set.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateJacobian()
    {
      using ModelT = BranchBreakers<ScalarT, IdxT>;
      using Fn     = Enzyme::Sparse::MemberFunctions;

      const auto n_y       = static_cast<size_t>(size_);
      const auto n_ws      = static_cast<size_t>(ws_.getSize());
      const auto n_out     = output_indices_.size();
      const auto bus1_rows = static_cast<size_t>(bus1_->size());
      const auto bus2_rows = static_cast<size_t>(bus2_->size());
      const auto bus1_cols = static_cast<size_t>(bus1_->y().getSize());
      const auto bus2_cols = static_cast<size_t>(bus2_->y().getSize());

      if (J_rows_buffer_ == nullptr)
      {
        const auto bus_size    = bus1_rows + bus2_rows;
        const auto buffer_size = 2 * n_y * n_y + n_y * n_ws
                                 + bus_size * (bus_size + n_y)
                                 + n_out * (n_y + wb_indices_.size());
        J_rows_buffer_ = new IdxT[buffer_size];
        J_cols_buffer_ = new IdxT[buffer_size];
        J_vals_buffer_ = new RealT[buffer_size];
      }

      // Bus-voltage columns in wb_ order; a bus without variables leaves INVALID columns
      std::ranges::fill(wb_indices_, INVALID_INDEX<IdxT>);
      std::ranges::copy(bus1_->getVariableIndices(), wb_indices_.begin());
      std::ranges::copy(bus2_->getVariableIndices(), wb_indices_.begin() + 2);

      ScalarT*    y    = y_.getData();
      ScalarT*    yp   = yp_.getData();
      ScalarT*    wb   = wb_.getData();
      ScalarT*    ws   = ws_.getData();
      const IdxT* res1 = bus1_->getResidualIndices().data();
      const IdxT* res2 = bus2_->getResidualIndices().data();
      const IdxT* var1 = bus1_->getVariableIndices().data();
      const IdxT* var2 = bus2_->getVariableIndices().data();

      nnz_ = 0;

      // Breaker latches
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

      // Bus currents: bus-k rows, bus-l voltage columns, as in Branch
      Enzyme::Sparse::DhDwb<ModelT, Fn::BusResidual11>::eval(this,
                                                             bus1_rows,
                                                             bus1_cols,
                                                             res1,
                                                             var1,
                                                             y,
                                                             yp,
                                                             wb,
                                                             J_rows_buffer_,
                                                             J_cols_buffer_,
                                                             J_vals_buffer_,
                                                             nnz_);

      Enzyme::Sparse::DhDwb<ModelT, Fn::BusResidual12>::eval(this,
                                                             bus1_rows,
                                                             bus2_cols,
                                                             res1,
                                                             var2,
                                                             y,
                                                             yp,
                                                             wb + 2,
                                                             J_rows_buffer_,
                                                             J_cols_buffer_,
                                                             J_vals_buffer_,
                                                             nnz_);

      Enzyme::Sparse::DhDwb<ModelT, Fn::BusResidual21>::eval(this,
                                                             bus2_rows,
                                                             bus1_cols,
                                                             res2,
                                                             var1,
                                                             y,
                                                             yp,
                                                             wb,
                                                             J_rows_buffer_,
                                                             J_cols_buffer_,
                                                             J_vals_buffer_,
                                                             nnz_);

      Enzyme::Sparse::DhDwb<ModelT, Fn::BusResidual22>::eval(this,
                                                             bus2_rows,
                                                             bus2_cols,
                                                             res2,
                                                             var2,
                                                             y,
                                                             yp,
                                                             wb + 2,
                                                             J_rows_buffer_,
                                                             J_cols_buffer_,
                                                             J_vals_buffer_,
                                                             nnz_);

      // ... and breaker-state columns
      Enzyme::Sparse::DhDy<ModelT, Fn::BusResidual11>::eval(this,
                                                            bus1_rows,
                                                            n_y,
                                                            res1,
                                                            variable_indices_.data(),
                                                            y,
                                                            yp,
                                                            wb,
                                                            J_rows_buffer_,
                                                            J_cols_buffer_,
                                                            J_vals_buffer_,
                                                            nnz_);

      Enzyme::Sparse::DhDy<ModelT, Fn::BusResidual12>::eval(this,
                                                            bus1_rows,
                                                            n_y,
                                                            res1,
                                                            variable_indices_.data(),
                                                            y,
                                                            yp,
                                                            wb + 2,
                                                            J_rows_buffer_,
                                                            J_cols_buffer_,
                                                            J_vals_buffer_,
                                                            nnz_);

      Enzyme::Sparse::DhDy<ModelT, Fn::BusResidual21>::eval(this,
                                                            bus2_rows,
                                                            n_y,
                                                            res2,
                                                            variable_indices_.data(),
                                                            y,
                                                            yp,
                                                            wb,
                                                            J_rows_buffer_,
                                                            J_cols_buffer_,
                                                            J_vals_buffer_,
                                                            nnz_);

      Enzyme::Sparse::DhDy<ModelT, Fn::BusResidual22>::eval(this,
                                                            bus2_rows,
                                                            n_y,
                                                            res2,
                                                            variable_indices_.data(),
                                                            y,
                                                            yp,
                                                            wb + 2,
                                                            J_rows_buffer_,
                                                            J_cols_buffer_,
                                                            J_vals_buffer_,
                                                            nnz_);

      // Computed-output gradients; SystemModel chains them into the readers' rows
      Enzyme::Sparse::DfDy<ModelT, Fn::Outputs>::eval(this,
                                                      n_out,
                                                      n_y,
                                                      output_indices_.data(),
                                                      variable_indices_.data(),
                                                      y,
                                                      yp,
                                                      wb,
                                                      ws,
                                                      J_rows_buffer_,
                                                      J_cols_buffer_,
                                                      J_vals_buffer_,
                                                      nnz_);

      Enzyme::Sparse::DfDwb<ModelT, Fn::Outputs>::eval(this,
                                                       n_out,
                                                       wb_indices_.size(),
                                                       output_indices_.data(),
                                                       wb_indices_.data(),
                                                       y,
                                                       yp,
                                                       wb,
                                                       ws,
                                                       J_rows_buffer_,
                                                       J_cols_buffer_,
                                                       J_vals_buffer_,
                                                       nnz_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class BranchBreakers<double, long int>;
    template class BranchBreakers<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
