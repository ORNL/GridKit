
#pragma once

namespace GridKit
{
  template <typename scalar_type>
  class ScalarTraits
  {
  };

  template <>
  class ScalarTraits<double>
  {
  public:
    typedef double RealT;
    typedef double NormT;
  };

} // namespace GridKit
