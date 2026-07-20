#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief By default, Jacobians are not available
     *
     * DependencyTracking::Variable stores the Jacobian as dependency maps. 
     * @todo Construct a Jacobian based on the dependency maps.
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      return false;
    }

    // Available template instantiations
    // template class SystemModel<DependencyTracking::Variable, long int>;
    template class SystemModel<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
