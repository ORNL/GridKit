/**
 * @file GenClassicalEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "GenClassicalImpl.hpp"

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
    int GenClassical<ScalarT, IdxT>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for GenClassical..." << std::endl;
      Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

      J_.zeroMatrix();

      GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::GenClassical<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual,
                                    ScalarT,
                                    IdxT>::eval(this,
                                                f_.size(),
                                                y_.size(),
                                                this->getResidualIndices(),
                                                this->getVariableIndices(),
                                                y_.data(),
                                                yp_.data(),
                                                wb_.data(),
                                                alpha_,
                                                J_);

      GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::GenClassical<ScalarT, IdxT>,
                                     GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual,
                                     ScalarT,
                                     IdxT>::eval(this,
                                                 f_.size(),
                                                 static_cast<size_t>(bus_->size()),
                                                 this->getResidualIndices(),
                                                 bus_->getVariableIndices(),
                                                 y_.data(),
                                                 yp_.data(),
                                                 (bus_->y()).data(),
                                                 J_);

      GridKit::Enzyme::Sparse::DhDy<GridKit::PhasorDynamics::GenClassical<ScalarT, IdxT>,
                                    GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                    ScalarT,
                                    IdxT>::eval(this,
                                                static_cast<size_t>(bus_->size()),
                                                y_.size(),
                                                bus_->getResidualIndices(),
                                                this->getVariableIndices(),
                                                y_.data(),
                                                yp_.data(),
                                                wb_.data(),
                                                J_);

      return 0;
    }

    // Available template instantiations
    template class GenClassical<double, long int>;
    template class GenClassical<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
