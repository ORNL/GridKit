/**
 * @file ComponentData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Data of EMT components.
 */

#pragma once

#include <array>
#include <tuple>
#include <vector>

#include <GridKit/Constants.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Parameter `key` of the model data, and its member
    template <typename KeyT, typename DataT, typename MemberT>
    struct Parameter
    {
      KeyT    key;
      MemberT DataT::* member;
    };

    /// Clang 16 has no aggregate deduction
    template <typename KeyT, typename DataT, typename MemberT>
    Parameter(KeyT, MemberT DataT::*) -> Parameter<KeyT, DataT, MemberT>;

    template <typename T>
    using Matrix = std::vector<std::vector<T>>;

    /// Real and imaginary parts
    template <typename RealT>
    using Pair = std::array<RealT, 2>;
  } // namespace EMT
} // namespace GridKit
