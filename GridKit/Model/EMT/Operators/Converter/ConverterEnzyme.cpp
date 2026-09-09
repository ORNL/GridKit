#include <GridKit/AutomaticDifferentiation/Enzyme/OutputGradient.hpp>

#include "ConverterImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::appendOutputGradient(
        Outputs output, typename SignalT::GradientT& gradient, RealT scale) const
    {
      if (verify() != 0 || static_cast<size_t>(output) >= output_port_.size())
        throw std::logic_error("Cannot differentiate an invalid Converter output");
      std::array<ScalarT, 4> values{};
      for (size_t n = 0; n < values.size(); ++n)
        values[n] = input_[n]->read();
      const auto partials = Enzyme::OutputGradient<Converter<ScalarT, IdxT>, 4>::eval(this, output, values);
      for (size_t n = 0; n < input_.size(); ++n)
        input_[n]->appendGradient(gradient, scale * partials[n]);
    }

    template class Converter<double, long int>;
    template class Converter<double, size_t>;
  } // namespace EMT
} // namespace GridKit
