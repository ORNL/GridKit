/**
 * @file GenrouEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

#include "GenrouImpl.hpp"

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
    int Genrou<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Genrou..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      GridKit::Enzyme::Sparse::InternalJacobianWithSignal<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                                          GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal,
                                                          ScalarT,
                                                          IdxT>::eval(this,
                                                                      f_.size(),
                                                                      y_.size(),
                                                                      this->getResidualIndices(),
                                                                      this->getVariableIndices(),
                                                                      y_.data(),
                                                                      yp_.data(),
                                                                      wb_.data(),
                                                                      ws_.data(),
                                                                      J_);

      J_.printMatrix("Genrou internal Jacobian");

      GridKit::Enzyme::Sparse::ExternalJacobian<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                                GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal,
                                                ScalarT,
                                                IdxT>::eval(this,
                                                            f_.size(),
                                                            ws_.size(),
                                                            this->getResidualIndices(),
                                                            ws_indices_,
                                                            y_.data(),
                                                            yp_.data(),
                                                            wb_.data(),
                                                            ws_.data(),
                                                            J_);

      J_.printMatrix("Genrou Jacobian after signal evaluation");

      GridKit::Enzyme::Sparse::BusJacobian<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                           GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                           ScalarT,
                                           IdxT>::eval(this,
                                                       static_cast<size_t>(bus_->size()),
                                                       static_cast<size_t>(bus_->size()),
                                                       bus_->getResidualIndices(),
                                                       bus_->getVariableIndices(),
                                                       y_.data(),
                                                       yp_.data(),
                                                       (bus_->y()).data(),
                                                       bus_->getJacobian());

      return 0;
    }

    // Available template instantiations
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
