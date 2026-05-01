/**
 * @file ArgVector.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 */

#pragma once

#include <array>
#include <string>
#include <vector>

#include <GridKit/Utilities/CliArgs/ArgValue.hpp>

namespace GridKit
{
  namespace Utilities
  {

    // Forward declaration
    struct CliArgsImpl;

    /**
     * @brief Represents a set of 0 or more values associated with a given
     * command-line option
     *
     * This class allows the internal data to be interpreted as a discrete set
     * of values (std::array) of a given type or as a single value, hopefully
     * making it more intuitive for users. For example...
     *
     * If an ArgVector v is expected to hold one value of type `double`, that
     * value can be extracted with `auto r = v.as<double>();`. This means `r`
     * will be a `double` value assigned with the contents of the first element
     * of `v`.
     *
     * If an ArgVector v is expected to hold two values of type `int`, on the
     * other hand, the values can be extracted as a `std::array<int, 2>` using
     * `as<int, 2>()`. Using `std::array` allows structured bindings for the
     * user: `auto [a, b] = v.as<int, 2>();`
     */
    class ArgVector
    {
    public:
      /**
       * @brief Default is empty (no values)
       */
      ArgVector() = default;

      /**
       * @brief Construct with a single value
       */
      template <typename T>
      ArgVector(T&& val)
        : vec{val}
      {
      }

      ArgVector(std::initializer_list<ArgValue> vals)
        : vec{vals}
      {
      }

      /**
       * @brief Interpret as `N` values of type `T`
       */
      template <typename T, std::size_t N>
      std::array<T, N> as() const
      {
        assert(vec.size() == N);
        std::array<T, N> ret;
        for (std::size_t i = 0; i < N; ++i)
        {
          ret[i] = vec[i].as<T>();
        }
        return ret;
      }

      /**
       * @brief Interpret as `N` values of type `std::string`
       */
      template <std::size_t N>
      decltype(auto) as() const
      {
        return as<std::string, N>();
      }

      /**
       * @brief Interpret as single value of type `T`
       */
      template <typename T>
      decltype(auto) as() const
      {
        return vec[0].as<T>();
      }

      /**
       * @brief Interpret as single `std::string` value
       */
      const std::string& operator()() const
      {
        return vec[0].get();
      }

      /**
       * @brief Get one of the contained `ArgValue` objects
       */
      const ArgValue& operator[](std::size_t i) const
      {
        return vec[i];
      }

      /**
       * @brief Check for existence of values
       */
      bool empty() const
      {
        return vec.empty();
      }

      /**
       * @brief Get number of `ArgValue` objects
       */
      std::size_t size() const
      {
        return vec.size();
      }

    private:
      /// Internal set of values
      std::vector<ArgValue> vec;

      friend struct CliArgsImpl;
    };

  } // namespace Utilities
} // namespace GridKit
