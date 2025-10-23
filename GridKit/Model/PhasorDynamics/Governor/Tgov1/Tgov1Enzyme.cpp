/**
 * @file Tgov1Enzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

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
        std::cout << "Evaluate Jacobian for Tgov1..." << std::endl;
        std::cout << "Jacobian evaluation is experimental!" << std::endl;

        GridKit::Enzyme::Sparse::InternalJacobian<GridKit::PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>,
                                                  GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual,
                                                  ScalarT,
                                                  IdxT>::eval(this,
                                                              f_.size(),
                                                              y_.size(),
                                                              this->getResidualIndices(),
                                                              this->getVariableIndices(),
                                                              y_.data(),
                                                              yp_.data(),
                                                              w_.data(),
                                                              J_);

        J_.printMatrix("Tgov1 internal Jacobian");

        return 0;
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
