/**
 * @file ArgValue.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 */

#pragma once

#include <sstream>
#include <string>

namespace GridKit
{
  namespace Utilities
  {

    /**
     * @brief Represents a single value of arbitrary type
     *
     * The value is internally represented as a std::string that is parsed when
     * requested as a specific type
     */
    class ArgValue
    {
      /**
       * @brief Check that T != ArgValue
       *
       */
      template <typename T>
      static constexpr bool notAnArgValue =
          !std::is_same_v<std::remove_reference_t<T>, ArgValue>;

      /**
       * @brief Check that T != std::initializer_list<ArgValue>
       */
      template <typename T>
      static constexpr bool notAnArgValueList =
          !std::is_same_v<std::remove_reference_t<T>, std::initializer_list<ArgValue>>;

      /**
       * @brief Check that T is not an ArgValue or list of ArgValues
       */
      template <typename T>
      static constexpr bool notArgValue = notAnArgValue<T> && notAnArgValueList<T>;

    public:
      /**
       * @brief Default construction results in empty value
       *
       * @note the SFINAE parameter is used to work around an issue with libc++
       */
      ArgValue() = default;

      /**
       * @brief Copy constructor
       */
      ArgValue(const ArgValue&) = default;

      /**
       * @brief Move constructor
       */
      ArgValue(ArgValue&&) = default;

      /**
       * @brief Construct from any value
       */
      template <typename T>
      ArgValue(T&& val, std::enable_if_t<notArgValue<T>, int> = 0)
        : value_((std::stringstream() << val).str())
      {
      }

      ArgValue& operator=(const ArgValue&) = default;
      ArgValue& operator=(ArgValue&&)      = default;

      /**
       * @brief Assign any value
       */
      template <typename T>
      ArgValue& operator=(T&& val)
      {
        value_ = (std::stringstream() << val).str();
        return *this;
      }

      /**
       * @brief Check if no value is contained
       */
      bool empty() const
      {
        return value_.empty();
      }

      /**
       * @brief Get string representation
       */
      const std::string& get() const
      {
        return value_;
      }

      /**
       * @brief Get string representation
       */
      const std::string& operator()() const
      {
        return get();
      }

      /**
       * @brief Get value of specific expected type
       */
      template <typename T>
      T as() const
      {
        T ret;
        std::stringstream(value_) >> ret;
        return ret;
      }

    private:
      /// Internal representation
      std::string value_{};
    };

  } // namespace Utilities
} // namespace GridKit
