#include "ConverterImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    void Converter<scalar_type, index_type>::appendOutputGradient(
        Outputs, typename SignalT::GradientT&, RealT) const
    {
      throw std::logic_error("Computed-output derivatives require Enzyme");
    }

    template class Converter<double, long int>;
    template class Converter<double, size_t>;
  } // namespace EMT
} // namespace GridKit
