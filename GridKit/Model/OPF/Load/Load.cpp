#include <cstddef>

#include "LoadImpl.hpp"

namespace GridKit
{
  namespace OPF
  {
    template class Load<double, int>;
    template class Load<double, long int>;
    template class Load<double, std::size_t>;
  } // namespace OPF
} // namespace GridKit
