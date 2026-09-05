#include "DependentVoltageSourceImpl.hpp"

namespace GridKit::EMT
{
  template class DependentVoltageSource<double, long int>;
  template class DependentVoltageSource<double, size_t>;
} // namespace GridKit::EMT
