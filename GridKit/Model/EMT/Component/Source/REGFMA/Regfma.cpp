#include "RegfmaImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Regfma<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      throw std::logic_error("Regfma requires an Enzyme-enabled build for Jacobian evaluation");
    }

    template class Regfma<double, long int>;
    template class Regfma<double, size_t>;
  } // namespace EMT
} // namespace GridKit
