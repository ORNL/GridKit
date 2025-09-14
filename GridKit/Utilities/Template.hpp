#pragma once

namespace GridKit
{
  namespace Utility
  {
    /**
     * @brief A simple struct for creating an overload visitor, such as for `std::visit`.
     * @see https://en.cppreference.com/w/cpp/utility/variant/visit2.html
     */
    template <class... Ts>
    struct OverloadVisitor : Ts...
    {
      using Ts::operator()...;
    };
  }; // namespace Utility
}; // namespace GridKit
