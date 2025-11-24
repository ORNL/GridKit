/**
 * @file BranchEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

#include "BranchImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Jacobian evaluation experimental
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateJacobian()
    {
      //std::cout << "Evaluate Jacobian for Branch..." << std::endl;
      //std::cout << "Jacobian evaluation is experimental!" << std::endl;

      J_.zeroMatrix();

      // Bus 1 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::BusJacobian<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                           GridKit::Enzyme::Sparse::MemberFunctions::BusResidual11,
                                           ScalarT,
                                           IdxT>::eval(this,
                                                       static_cast<size_t>(bus1_->size()),
                                                       static_cast<size_t>((bus1_->y()).size()),
                                                       bus1_->getResidualIndices(),
                                                       bus1_->getVariableIndices(),
                                                       y_.data(),
                                                       yp_.data(),
                                                       (bus1_->y()).data(),
                                                       bus1_->getJacobian());

      // Bus 2 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::BusJacobian<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                           GridKit::Enzyme::Sparse::MemberFunctions::BusResidual22,
                                           ScalarT,
                                           IdxT>::eval(this,
                                                       static_cast<size_t>(bus2_->size()),
                                                       static_cast<size_t>((bus2_->y()).size()),
                                                       bus2_->getResidualIndices(),
                                                       bus2_->getVariableIndices(),
                                                       y_.data(),
                                                       yp_.data(),
                                                       (bus2_->y()).data(),
                                                       bus2_->getJacobian());

      // Off-diagonal Jacobian block (Bus2 variables) owned by the branch
      GridKit::Enzyme::Sparse::BranchJacobian<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                              GridKit::Enzyme::Sparse::MemberFunctions::BusResidual12,
                                              ScalarT,
                                              IdxT>::eval(this,
                                                          static_cast<size_t>(bus1_->size()),
                                                          static_cast<size_t>((bus2_->y()).size()),
                                                          bus1_->getResidualIndices(),
                                                          bus2_->getVariableIndices(),
                                                          y_.data(),
                                                          yp_.data(),
                                                          (bus2_->y()).data(),
                                                          J_);

      // Off-diagonal Jacobian block (Bus1 variables) owned by the branch
      GridKit::Enzyme::Sparse::BranchJacobian<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                              GridKit::Enzyme::Sparse::MemberFunctions::BusResidual21,
                                              ScalarT,
                                              IdxT>::eval(this,
                                                          static_cast<size_t>(bus2_->size()),
                                                          static_cast<size_t>((bus1_->y()).size()),
                                                          bus2_->getResidualIndices(),
                                                          bus1_->getVariableIndices(),
                                                          y_.data(),
                                                          yp_.data(),
                                                          (bus1_->y()).data(),
                                                          J_);

      return 0;
    }

    // Available template instantiations
    template class Branch<double, long int>;
    template class Branch<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
