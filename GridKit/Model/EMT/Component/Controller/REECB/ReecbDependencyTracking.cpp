
#include "ReecbImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      /**
       * @brief Jacobian evaluation not implemented
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::assembleJacobian(RealT, RealT)
      {
        Log::misc() << "Evaluate Jacobian for Reecb..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

        return 0;
      }

      // Available template instantiations
      template class Reecb<DependencyTracking::Variable, long int>;
      template class Reecb<DependencyTracking::Variable, size_t>;

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
