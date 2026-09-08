
#include "OuterPowerControlImpl.hpp"

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
      int OuterPowerControl<scalar_type, index_type>::assembleJacobian(RealT, RealT)
      {
        Log::misc() << "Evaluate Jacobian for OuterPowerControl..." << std::endl;
        Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

        return 0;
      }

      // Available template instantiations
      template class OuterPowerControl<double, long int>;
      template class OuterPowerControl<double, size_t>;

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
