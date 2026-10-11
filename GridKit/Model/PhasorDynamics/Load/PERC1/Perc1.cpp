#include "Perc1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Report that a separate Jacobian is unavailable in the plain build.
     */
    template <typename scalar_type, typename index_type>
    int Perc1<scalar_type, index_type>::evaluateJacobian()
    {
      Log::misc() << "Perc1: Jacobian evaluation is not implemented\n";
      return 0;
    }

    template class Perc1<double, long int>;
    template class Perc1<double, size_t>;
  } // namespace PhasorDynamics
} // namespace GridKit
