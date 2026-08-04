#include <cstddef>

#include "BusImpl.hpp"

namespace GridKit
{
  namespace OPF
  {
    template class Bus<double, int>;
    template class Bus<double, long int>;
    template class Bus<double, std::size_t>;
  } // namespace OPF
} // namespace GridKit
