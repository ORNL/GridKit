#include <GridKit/AutomaticDifferentiation/Enzyme/OutputGradient.hpp>

#include "PwmImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      void Pwm<scalar_type, index_type>::appendOutputGradient(
          Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
      {
        if (verify() != 0 || static_cast<size_t>(output) >= output_port_.size())
          throw std::logic_error("Cannot differentiate an invalid Pwm output");
        if (!hasInput())
          return;
        const auto values   = inputValues();
        const auto partials = Enzyme::OutputGradient<Pwm<ScalarT, IdxT>, 4>::eval(this, output, values);
        for (size_t n = 0; n < input_.size(); ++n)
          input_[n]->appendGradient(gradient, scale * partials[n]);
      }

      template class Pwm<double, long int>;
      template class Pwm<double, size_t>;
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
