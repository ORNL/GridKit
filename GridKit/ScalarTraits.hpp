
#pragma once

namespace GridKit
{
  template <class ScalarT>
  class ScalarTraits
  {
  };

  template <>
  class ScalarTraits<double>
  {
  public:
    typedef double real_type;
    typedef double norm_type;
    typedef double scalar_type;
  };

} // namespace GridKit
