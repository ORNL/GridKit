#include "SystemModelImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief By default, Jacobians are not available
     *
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      Log::warning() << "GridKit was not built with Enzyme. "
                     << "Falling back to dense Jacobians in PhasorDynamics.\n";

      return false;
    }

    // Available template instantiations
    // template class SystemModel<double, long int>;
    template class SystemModel<double, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit
