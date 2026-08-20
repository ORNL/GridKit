#pragma once

#include <ranges>

#include <magic_enum/magic_enum.hpp>

namespace GridKit
{
  namespace Utilities
  {

    /// Helper for checking contiguity of an enum.
    template <typename T>
    consteval auto is_contiguous_enum() -> bool
    {
      if constexpr (magic_enum::enum_count<T>() == 0)
      {
        return true;
      }
      else
      {
        std::size_t index = 0;
        for (auto variant : magic_enum::enum_values<T>())
        {
          if (static_cast<std::size_t>(variant) != index++)
          {
            return false;
          }
        }

        return true;
      }
    };

    /// An enum with a contiguous sequence of variants, a known maximum value,
    /// and an underlying type of @ref std::size_t.
    template <typename T>
    concept SizedEnum = std::is_enum_v<T>
                        && std::is_same_v<std::underlying_type_t<T>,
                                          std::size_t>
                        && is_contiguous_enum<T>();

    /// Simplest enum satisfying [`SizedEnum`].
    enum class EmptyEnum : std::size_t
    {
    };

    /// Variant of `magic_enum::enum_cast` that permits passing empty enums.
    template <SizedEnum T>
    constexpr auto enum_cast(std::string_view value) -> std::optional<T>
    {
      if constexpr (magic_enum::enum_count<T>() == 0)
      {
        return std::nullopt;
      }
      else
      {
        return magic_enum::enum_cast<T>(value);
      }
    }

    /// Variant of `magic_enum::enum_cast` that permits passing empty enums.
    template <SizedEnum T, typename BinaryPredicate = std::equal_to<>>
    constexpr auto enum_cast(std::string_view value, BinaryPredicate predicate)
        -> std::optional<T>
    {
      if constexpr (magic_enum::enum_count<T>() == 0)
      {
        return std::nullopt;
      }
      else
      {
        return magic_enum::enum_cast<T>(value, predicate);
      }
    }

    /// Upper bound on attainable values of the specified @ref SizedEnum.
    template <SizedEnum T>
    consteval auto enum_size() -> std::size_t
    {
      return magic_enum::enum_count<T>();
    }

    /// Check if a @ref SizedEnum value corresponds to a variant.
    template <SizedEnum T>
    constexpr auto contained_within(T v) -> bool
    {
      return static_cast<std::size_t>(v) < enum_size<T>();
    }

    /// Enumerate over the variants of a @ref SizedEnum.
    template <SizedEnum T>
    constexpr auto enum_values()
    {
      return std::views::iota(
                 std::size_t{0},
                 static_cast<std::size_t>(magic_enum::enum_count<T>()))
             | std::views::transform([](std::size_t i)
                                     { return static_cast<T>(i); });
    }
  } // namespace Utilities
} // namespace GridKit
