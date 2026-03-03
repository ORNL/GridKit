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
    public:
      /**
       * @brief Default construction results in empty value
       */
      ArgValue() = default;

      /**
       * @brief Construct from any value
       */
      template <typename T>
      ArgValue(T&& val)
        : value_((std::stringstream() << val).str())
      {
      }

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
