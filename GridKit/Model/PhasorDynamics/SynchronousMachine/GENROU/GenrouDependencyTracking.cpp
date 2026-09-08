#include "GenrouImpl.hpp"

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
    int Genrou<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate DependencyTracking Jacobian for Genrou...\n";
      Log::misc() << "Jacobian evaluation is experimental!\n";

      this->constructCsr();

      return 0;
    }

    // Available template instantiations
    template class Genrou<DependencyTracking::Variable, long int>;
    template class Genrou<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
