#pragma once

namespace GridKit
{
  namespace Utilities
  {
    template <typename TEnum>
    concept SizedEnum = requires {
      std::is_enum_v<TEnum>;
      TEnum::SIZE;
    };

    template <SizedEnum TEnum>
    inline constexpr std::underlying_type_t<TEnum> enumId(TEnum v) noexcept
    {
      return static_cast<std::underlying_type_t<TEnum>>(v);
    }

    template <SizedEnum TEnum>
    inline constexpr auto enumSize = enumId(TEnum::SIZE);
  } // namespace Utilities
} // namespace GridKit
