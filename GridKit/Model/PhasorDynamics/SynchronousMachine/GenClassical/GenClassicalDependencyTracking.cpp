
#include "GenClassicalImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Evaluate DependencyTracking::Variable Jacobian.
     *
     * @note Currently only used for testing.
     *
     * DependencyTracking::Variable stores the Jacobian as dependency maps,
     * updated during calls to evaluateResidual().
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate DependencyTracking Jacobian for GenClassical...\n";
      Log::misc() << "Jacobian evaluation is experimental!\n";

      this->constructCsr();

      return 0;
    }

    // Available template instantiations
    template class GenClassical<DependencyTracking::Variable, long int>;
    template class GenClassical<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
