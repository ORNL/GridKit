
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
    typedef double RealT;
    typedef double NormT;
  };

} // namespace GridKit
