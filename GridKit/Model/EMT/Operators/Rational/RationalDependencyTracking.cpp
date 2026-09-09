/**
 * @file RationalDependencyTracking.cpp
 * @brief Template instantiations for the EMT Rational model.
 */

#include "RationalImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::assembleJacobian(RealT, RealT)
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::appendOutputGradient(
        IdxT, typename SignalT::GradientT&, RealT) const
    {
      throw std::logic_error("Rational output gradients require Enzyme");
    }

    template class Rational<DependencyTracking::Variable, long int>;
    template class Rational<DependencyTracking::Variable, size_t>;
  } // namespace EMT
} // namespace GridKit
