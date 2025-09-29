
#include "LoadImpl.hpp"
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
    int Load<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Load..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      GridKit::Enzyme::Sparse::BusJacobian<GridKit::PhasorDynamics::Load<ScalarT, IdxT>,
                                           GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                           ScalarT,
                                           IdxT>::eval(
          this,
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
    template class Load<double, long int>;
    template class Load<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
