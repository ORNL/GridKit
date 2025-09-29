
#include "GenrouImpl.hpp"
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
    int Genrou<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Genrou..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      GridKit::Enzyme::Sparse::InternalJacobian<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
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

      J_.printMatrix("Genrou internal Jacobian");

      GridKit::Enzyme::Sparse::BusJacobian<GridKit::PhasorDynamics::Genrou<ScalarT, IdxT>,
                                           GridKit::Enzyme::Sparse::MemberFunctions::BusResidual,
                                           ScalarT,
                                           IdxT>(
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
    template class Genrou<double, long int>;
    template class Genrou<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
