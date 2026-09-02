
#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Governor
    {
      /**
       * @brief Jacobian evaluation not implemented
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Tgov1..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

        return 0;
      }

      // Available template instantiations
      template class Tgov1<DependencyTracking::Variable, long int>;
      template class Tgov1<DependencyTracking::Variable, size_t>;

    } // namespace Governor
  } // namespace EMT
} // namespace GridKit
