#include "PwmImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::appendOutputGradient(
          Outputs, typename SignalT::GradientT&, RealT) const
      {
        throw std::logic_error("Computed-output derivatives require Enzyme");
      }

      template class Pwm<double, long int>;
      template class Pwm<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
