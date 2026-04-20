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
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateJacobian()
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
                                                 static_cast<size_t>(this->port(0)->size()),
                                                 static_cast<size_t>((this->port(0)->y()).size()),
                                                 (this->port(0)->getResidualIndices()).data(),
                                                 (this->port(0)->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (this->port(0)->y()).data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 this->port(0)->getJacobian());

      // Bus 2 diagonal Jacobian block owned by the bus
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual22,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(this->port(1)->size()),
                                                 static_cast<size_t>((this->port(1)->y()).size()),
                                                 (this->port(1)->getResidualIndices()).data(),
                                                 (this->port(1)->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (this->port(1)->y()).data(),
                                                 J_rows_buffer_,
                                                 J_cols_buffer_,
                                                 J_vals_buffer_,
                                                 this->port(1)->getJacobian());

      // Off-diagonal Jacobian block (Bus2 variables) owned by the branch
      GridKit::Enzyme::Sparse::DhDwb<GridKit::PhasorDynamics::Branch<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::BusResidual12,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 static_cast<size_t>(this->port(0)->size()),
                                                 static_cast<size_t>((this->port(1)->y()).size()),
                                                 (this->port(0)->getResidualIndices()).data(),
                                                 (this->port(1)->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (this->port(1)->y()).data(),
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
                                                 static_cast<size_t>(this->port(1)->size()),
                                                 static_cast<size_t>((this->port(0)->y()).size()),
                                                 (this->port(1)->getResidualIndices()).data(),
                                                 (this->port(0)->getVariableIndices()).data(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (this->port(0)->y()).data(),
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
