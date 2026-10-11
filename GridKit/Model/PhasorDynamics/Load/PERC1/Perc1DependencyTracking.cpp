#include "Perc1Impl.hpp"

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
     */
    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::evaluateJacobian()
    {
      this->constructCsr();

      return 0;
    }

    // Available template instantiations
    template class Perc1<DependencyTracking::Variable, long int>;
    template class Perc1<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
