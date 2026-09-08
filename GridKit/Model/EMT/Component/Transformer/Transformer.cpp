
#include "TransformerImpl.hpp"

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
    int Transformer<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      Log::misc() << "Evaluate Jacobian for Transformer..." << std::endl;
      Log::misc() << "Jacobian evaluation is not implemented!" << std::endl;

      return 0;
    }

    // Available template instantiations
    template class Transformer<double, long int>;
    template class Transformer<double, size_t>;

  } // namespace EMT
} // namespace GridKit
