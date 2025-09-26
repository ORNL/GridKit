
#include "BusFaultImpl.hpp"
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
    int BusFault<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for BusFault..." << std::endl;
      std::cout << "Jacobian evaluation is experimental!" << std::endl;

      if (status_)
      {
        GridKit::Enzyme::Sparse::BusJacobian<GridKit::PhasorDynamics::BusFault<ScalarT, IdxT>,
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
      }

      return 0;
    }

    // Available template instantiations
    template class BusFault<double, long int>;
    template class BusFault<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
