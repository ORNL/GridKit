/**
 * @file Tgov1Enzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Jacobian evaluation experimental
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int Tgov1<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Tgov1..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>,
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
                                                  alpha_,
                                                  J_);

        GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>,
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

        return 0;
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
