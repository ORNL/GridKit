// The affine residual has the same exact Jacobian with and without Enzyme.
#include "SignalDelayImpl.hpp"

template class GridKit::PhasorDynamics::SignalDelay<double, long int>;
template class GridKit::PhasorDynamics::SignalDelay<double, size_t>;
