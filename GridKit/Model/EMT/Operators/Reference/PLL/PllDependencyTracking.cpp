
#include "PllImpl.hpp"

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
    int Pll<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      Log::misc() << "Evaluate Jacobian for Pll..." << std::endl;
      Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class Pll<DependencyTracking::Variable, long int>;
    template class Pll<DependencyTracking::Variable, size_t>;

  } // namespace EMT
} // namespace GridKit
