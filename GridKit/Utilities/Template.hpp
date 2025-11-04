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

    template <class T>
    class Immovable
    {
    private:
      T inner_;

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

      T unwrap() &&
      {
        return std::move(inner_);
      }
    };
  }; // namespace Utility
}; // namespace GridKit
