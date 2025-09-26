
#include "GenClassicalImpl.hpp"
#include <AutomaticDifferentiation/Enzyme/SparseWrapper.hpp>

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

      GridKit::Enzyme::Sparse::ModelJacobian<GridKit::PhasorDynamics::GenClassical<ScalarT, IdxT>,
                                             GridKit::Enzyme::Sparse::MemberFunctions::InternalResidual,
                                             ScalarT,
                                             IdxT>(
          this,
          f_.size(),
          y_.size(),
          this->getResidualIndices(),
          this->getVariableIndices(),
          y_.data(),
          yp_.data(),
          w_.data(),
          J_);

      J_.printMatrix("GenClassical internal Jacobian");

      return 0;
    }

    // Available template instantiations
    template class GenClassical<double, long int>;
    template class GenClassical<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
