#pragma once

namespace GridKit
{
  namespace Utilities
  {
    template <typename T>
    concept SizedEnum = std::is_enum_v<T> && requires { T::SIZE; };

    enum class EmptyEnum : std::size_t
    {
      SIZE = 0
    };

    template <SizedEnum enum_type>
    constexpr std::size_t enum_size()
    {
      return static_cast<std::size_t>(enum_type::SIZE);
    }
  } // namespace Utilities
} // namespace GridKit
