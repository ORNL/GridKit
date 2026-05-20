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

      J_.zeroMatrix();
      if (J_rows_buffer_ == nullptr)
      {
        J_rows_buffer_ = new IdxT[4];
        J_cols_buffer_ = new IdxT[4];
        J_vals_buffer_ = new RealT[4];
      }

      // Bus 1 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual11,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(bus1_->size()),
                                                 static_cast<size_t>((bus1_->y()).size()),
                                                 (bus1_->getResidualIndices()).data(),
                                                 (bus1_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (bus1_->y()).data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 bus1_->getJacobian());

      // Bus 2 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual22,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(bus2_->size()),
                                                 static_cast<size_t>((bus2_->y()).size()),
                                                 (bus2_->getResidualIndices()).data(),
                                                 (bus2_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (bus2_->y()).data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 bus2_->getJacobian());

      // Off-diagonal Jacobian block (Bus2 variables) owned by the branch
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual12,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(bus1_->size()),
                                                 static_cast<size_t>((bus2_->y()).size()),
                                                 (bus1_->getResidualIndices()).data(),
                                                 (bus2_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (bus2_->y()).data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 J_,
                                                 true);

      // Off-diagonal Jacobian block (Bus1 variables) owned by the branch
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual21,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(bus2_->size()),
                                                 static_cast<size_t>((bus1_->y()).size()),
                                                 (bus2_->getResidualIndices()).data(),
                                                 (bus1_->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (bus1_->y()).data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 J_,
                                                 true);

      return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
