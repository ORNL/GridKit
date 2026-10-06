/**
 * @file BranchBreakersEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the breaker-terminated branch model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BranchBreakersImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Sparse Jacobian of the breaker-terminated branch
     *
     * Bus-1 then bus-2 index arrays address the two terminals; a size-0 bus
     * leaves INVALID_INDEX entries, which the sparse store skips.
     *
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateJacobian()
    {
      wb_indices_.assign(4, INVALID_INDEX<IdxT>);
      h_indices_.assign(4, INVALID_INDEX<IdxT>);
      std::ranges::copy(bus1_->getVariableIndices(), wb_indices_.begin());
      std::ranges::copy(bus2_->getVariableIndices(), wb_indices_.begin() + 2);
      std::ranges::copy(bus1_->getResidualIndices(), h_indices_.begin());
      std::ranges::copy(bus2_->getResidualIndices(), h_indices_.begin() + 2);

      const auto size        = static_cast<size_t>(size_);
      const auto bus_size    = wb_indices_.size();
      const auto signal_size = static_cast<size_t>(ws_.getSize());

      if (J_rows_buffer_ == nullptr)
      {
        auto buffer_size = 2 * size * size + size * signal_size + bus_size * size + bus_size * bus_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      using ModelT = GridKit::PhasorDynamics::BranchBreakers<scalar_type, index_type>;
      using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

      nnz_ = 0;

      GridKit::Enzyme::Sparse::DfDy<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                  size,
                                                                                  size,
                                                                                  (this->getResidualIndices()).data(),
                                                                                  (this->getVariableIndices()).data(),
                                                                                  y_.getData(),
                                                                                  yp_.getData(),
                                                                                  wb_.getData(),
                                                                                  ws_.getData(),
                                                                                  J_rows_buffer_,
                                                                                  J_cols_buffer_,
                                                                                  J_vals_buffer_,
                                                                                  nnz_);

      GridKit::Enzyme::Sparse::DfDyp<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                   size,
                                                                                   size,
                                                                                   (this->getResidualIndices()).data(),
                                                                                   (this->getVariableIndices()).data(),
                                                                                   y_.getData(),
                                                                                   yp_.getData(),
                                                                                   wb_.getData(),
                                                                                   ws_.getData(),
                                                                                   alpha_,
                                                                                   J_rows_buffer_,
                                                                                   J_cols_buffer_,
                                                                                   J_vals_buffer_,
                                                                                   nnz_);

      GridKit::Enzyme::Sparse::DfDws<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                   size,
                                                                                   signal_size,
                                                                                   (this->getResidualIndices()).data(),
                                                                                   ws_indices_.data(),
                                                                                   y_.getData(),
                                                                                   yp_.getData(),
                                                                                   wb_.getData(),
                                                                                   ws_.getData(),
                                                                                   J_rows_buffer_,
                                                                                   J_cols_buffer_,
                                                                                   J_vals_buffer_,
                                                                                   nnz_);

      GridKit::Enzyme::Sparse::DhDy<ModelT, Fn::BusResidual>::eval(this,
                                                                   bus_size,
                                                                   size,
                                                                   h_indices_.data(),
                                                                   (this->getVariableIndices()).data(),
                                                                   y_.getData(),
                                                                   yp_.getData(),
                                                                   wb_.getData(),
                                                                   J_rows_buffer_,
                                                                   J_cols_buffer_,
                                                                   J_vals_buffer_,
                                                                   nnz_);

      GridKit::Enzyme::Sparse::DhDwb<ModelT, Fn::BusResidual>::eval(this,
                                                                    bus_size,
                                                                    bus_size,
                                                                    h_indices_.data(),
                                                                    wb_indices_.data(),
                                                                    y_.getData(),
                                                                    yp_.getData(),
                                                                    wb_.getData(),
                                                                    J_rows_buffer_,
                                                                    J_cols_buffer_,
                                                                    J_vals_buffer_,
                                                                    nnz_);

      this->constructCoo();

      return 0;
    }

    template class BranchBreakers<double, long int>;
    template class BranchBreakers<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
