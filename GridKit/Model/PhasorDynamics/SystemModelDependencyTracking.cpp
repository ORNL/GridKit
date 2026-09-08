#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief By default, Jacobians are not available
     *
     * DependencyTracking::Variable stores the Jacobian as dependency maps,
     * updated during calls to evaluateResidual().
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      Log::warning() << "DependencyTracking::Variable Jacobians are only available for testing.\n"
                     << "Falling back to dense Jacobians for PhasorDyanmics simulations.\n";

      return false;
    }

    /**
     * @brief Evaluate system DependencyTracking::Variable Jacobian.
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      this->constructCsr();

      // Log::misc() << "System DependencyTracking Jacobian\n";
      // csr_jac_->print(Log::misc());

      return 0;
    }

    // Available template instantiations
    // template class SystemModel<DependencyTracking::Variable, long int>;
    template class SystemModel<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
