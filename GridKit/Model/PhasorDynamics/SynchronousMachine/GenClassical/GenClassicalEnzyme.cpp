/**
 * @file GenClassicalEnzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

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
      std::cout << "Evaluate Jacobian for GenClassical..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      GridKit::Enzyme::Sparse::InternalJacobian<GridKit::PhasorDynamics::GenClassical<ScalarT, IdxT>,
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
                                                            J_);

      J_.printMatrix("GenClassical internal Jacobian");

      return 0;
    }

    // Available template instantiations
    template class GenClassical<double, long int>;
    template class GenClassical<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
