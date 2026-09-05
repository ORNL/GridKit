#include "ConverterImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Converter<DependencyTracking::Variable, long int>;
    template class Converter<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
