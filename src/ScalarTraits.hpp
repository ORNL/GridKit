
#ifndef _SCALAR_TRAITS_HPP_
#define _SCALAR_TRAITS_HPP_

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

#endif // _SCALAR_TRAITS_HPP_
