
#include "VectorFitImpl.hpp"

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
    int VectorFit<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Evaluate Jacobian for VectorFit..." << std::endl;
      Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class VectorFit<DependencyTracking::Variable, long int>;
    template class VectorFit<DependencyTracking::Variable, size_t>;

  } // namespace EMT
} // namespace GridKit
