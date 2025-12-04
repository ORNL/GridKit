#pragma once

/**
 * @file
 */

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Requires `TEnum` to be an `enum` and have a `SIZE` member
     */
    template <typename TEnum>
    concept SizedEnum = requires {
      std::is_enum_v<TEnum>;
      TEnum::SIZE;
    };

    /**
     * @brief Convert an enum value to its underlying integer type
     */
    template <SizedEnum TEnum>
    inline constexpr std::underlying_type_t<TEnum> enumId(TEnum v) noexcept
    {
      return static_cast<std::underlying_type_t<TEnum>>(v);
    }

    /// Length of an enum
    template <SizedEnum TEnum>
    inline constexpr auto enumSize = enumId(TEnum::SIZE);
  } // namespace Utilities
} // namespace GridKit
