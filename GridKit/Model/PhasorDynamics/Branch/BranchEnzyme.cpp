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
        auto bus1_size   = static_cast<size_t>(bus1_->size());
        auto bus2_size   = static_cast<size_t>(bus2_->size());
        auto buffer_size = (bus1_size + bus2_size) * (bus1_size + bus2_size);
        J_rows_buffer_   = new IdxT[buffer_size];
        J_cols_buffer_   = new IdxT[buffer_size];
        J_vals_buffer_   = new RealT[buffer_size];
      }

      nnz_ = 0;

      // Bus 1 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual11>::eval(this,
                                                                                                    static_cast<size_t>(bus1_->size()),
                                                                                                    static_cast<size_t>((bus1_->y()).getSize()),
                                                                                                    (bus1_->getResidualIndices()).data(),
                                                                                                    (bus1_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus1_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      // Bus 2 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual22>::eval(this,
                                                                                                    static_cast<size_t>(bus2_->size()),
                                                                                                    static_cast<size_t>((bus2_->y()).getSize()),
                                                                                                    (bus2_->getResidualIndices()).data(),
                                                                                                    (bus2_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus2_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      // Off-diagonal Jacobian block (Bus2 variables) owned by the branch
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual12>::eval(this,
                                                                                                    static_cast<size_t>(bus1_->size()),
                                                                                                    static_cast<size_t>((bus2_->y()).getSize()),
                                                                                                    (bus1_->getResidualIndices()).data(),
                                                                                                    (bus2_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus2_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      // Off-diagonal Jacobian block (Bus1 variables) owned by the branch
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual21>::eval(this,
                                                                                                    static_cast<size_t>(bus2_->size()),
                                                                                                    static_cast<size_t>((bus1_->y()).getSize()),
                                                                                                    (bus2_->getResidualIndices()).data(),
                                                                                                    (bus1_->getVariableIndices()).data(),
                                                                                                    y_.getData(),
                                                                                                    yp_.getData(),
                                                                                                    bus1_->y().getData(),
                                                                                                    J_rows_buffer_,
                                                                                                    J_cols_buffer_,
                                                                                                    J_vals_buffer_,
                                                                                                    nnz_);

      this->constructCoo();

      return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
