
#include "LoadZImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Jacobian evaluation not implemented
     *
     * @return int - error code, 0 = success
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for LoadZ..." << std::endl;
      Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class LoadZ<DependencyTracking::Variable, long int>;
    template class LoadZ<DependencyTracking::Variable, size_t>;

  } // namespace EMT
} // namespace GridKit
