/**
 * @file Norton.cpp
 * @brief Template instantiations for the EMT Norton model.
 */

#include "NortonImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Norton<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      throw std::logic_error("Norton Jacobian requires Enzyme");
    }

    template class Norton<double, long int>;
    template class Norton<double, size_t>;
  } // namespace EMT
} // namespace GridKit
