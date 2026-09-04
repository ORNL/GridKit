/**
 * @file Signal model implementation.
 */
#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>

#include "SignalImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Signal<DependencyTracking::Variable, size_t>;
    template class Signal<DependencyTracking::Variable, long>;

  } // namespace EMT
} // namespace GridKit
