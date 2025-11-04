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

    /**
     * @brief A wrapper class which ensures that an underlying value is never moved.
     *         Useful if you'd like to take a reference to some value and you aren't sure whether
     *         or not that value can be moved later, invalidating that reference.
     * @tparam T Type of the wrapped value
     */
    template <class T>
    class Immovable
    {
    private:
      T inner_;

      // Simply delete move constructors
      Immovable(Immovable&& other)            = delete;
      Immovable& operator=(Immovable&& other) = delete;

    public:
      Immovable(T&& inner) : inner_(std::move(inner))
      {
      }

      operator T&()
      {
        return inner_;
      }

      operator T&() const
      {
        return inner_;
      }

      T& get()
      {
        return inner_;
      }

      const T& get() const
      {
        return inner_;
      }

      /**
       * @brief "Unwrap" the inner object, allowing for it to be moved again.
       */
      T unwrap() &&
      {
        return std::move(inner_);
      }
    };
  }; // namespace Utility
}; // namespace GridKit
