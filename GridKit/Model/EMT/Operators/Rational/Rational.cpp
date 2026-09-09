/**
 * @file Rational.cpp
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
      throw std::logic_error("Rational Jacobian requires Enzyme");
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::appendOutputGradient(
        IdxT, typename SignalT::GradientT&, RealT) const
    {
      throw std::logic_error("Rational output gradients require Enzyme");
    }

    template class Rational<double, long int>;
    template class Rational<double, size_t>;
  } // namespace EMT
} // namespace GridKit
