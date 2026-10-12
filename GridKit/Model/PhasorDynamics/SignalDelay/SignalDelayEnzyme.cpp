// The two affine Jacobian entries are exact and shared with the ordinary build.
#include "SignalDelayImpl.hpp"

template class GridKit::PhasorDynamics::SignalDelay<double, long int>;
template class GridKit::PhasorDynamics::SignalDelay<double, size_t>;
