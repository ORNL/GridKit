/**
 * @file BranchEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for Branch..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      if (J_rows_buffer_ == nullptr)
      {
        // Reserve space for the dense blocks.
        // The size of the buffer is the sum of maximum capacities of the blocks.
        // Enyme will compute the appropriate nnz from sparsification.
        auto size        = static_cast<size_t>(size_);
        auto bus1_size   = static_cast<size_t>(bus1_->size());
        auto bus2_size   = static_cast<size_t>(bus2_->size());
        auto buffer_size = size * (bus1_size + bus2_size) + size + bus1_size + bus2_size;
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      const auto ir1 = static_cast<size_t>(BranchInternalVariables::IR1);
      const auto ii1 = static_cast<size_t>(BranchInternalVariables::II1);
      const auto ir2 = static_cast<size_t>(BranchInternalVariables::IR2);
      const auto ii2 = static_cast<size_t>(BranchInternalVariables::II2);

      // Terminal-1 current equations with respect to bus-1 voltage variables
      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual11>::eval(this,
                                                                                                    2,
                                                                                                    static_cast<size_t>((bus1_->y()).getSize()),
                                                                                                    residual_indices_.data() + ir1,
                                                                                                    (bus1_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus1_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      // Terminal-2 current equations with respect to bus-2 voltage variables
      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual22>::eval(this,
                                                                                                    2,
                                                                                                    static_cast<size_t>((bus2_->y()).getSize()),
                                                                                                    residual_indices_.data() + ir2,
                                                                                                    (bus2_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus2_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      // Terminal-1 current equations with respect to bus-2 voltage variables
      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual12>::eval(this,
                                                                                                    2,
                                                                                                    static_cast<size_t>((bus2_->y()).getSize()),
                                                                                                    residual_indices_.data() + ir1,
                                                                                                    (bus2_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus2_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      // Terminal-2 current equations with respect to bus-1 voltage variables
      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual21>::eval(this,
                                                                                                    2,
                                                                                                    static_cast<size_t>((bus1_->y()).getSize()),
                                                                                                    residual_indices_.data() + ir2,
                                                                                                    (bus1_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus1_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      auto append_entry = [this](IdxT row, IdxT col, RealT value)
      {
        const auto index       = static_cast<size_t>(nnz_);
        J_rows_buffer_[index]  = row;
        J_cols_buffer_[index]  = col;
        J_vals_buffer_[index]  = value;
        nnz_                  += 1;
      };

      // Algebraic current equations: d(-I)/dI = -1.
      for (IdxT i = 0; i < size_; ++i)
      {
        append_entry(residual_indices_[static_cast<size_t>(i)],
                     variable_indices_[static_cast<size_t>(i)],
                     RealT{-1.0});
      }

      // Bus current-balance equations receive the branch terminal currents.
      if (bus1_->size() > 0)
      {
        append_entry(bus1_->getResidualIndices()[0], variable_indices_[ir1], RealT{1.0});
        append_entry(bus1_->getResidualIndices()[1], variable_indices_[ii1], RealT{1.0});
      }
      if (bus2_->size() > 0)
      {
        append_entry(bus2_->getResidualIndices()[0], variable_indices_[ir2], RealT{1.0});
        append_entry(bus2_->getResidualIndices()[1], variable_indices_[ii2], RealT{1.0});
      }

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
