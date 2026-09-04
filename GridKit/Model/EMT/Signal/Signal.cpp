/**
 * @file Signal model implementation.
 */
#include "SignalImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template class Signal<double, size_t>;
    template class Signal<double, long>;

  } // namespace EMT
} // namespace GridKit
