#pragma once

#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Convert a string to all uppercase
     */
    inline std::string toUpper(std::string str)
    {
      std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c)
                     { return std::toupper(c); });
      return str;
    }

    /**
     * @brief Attempt to parse a string to the given type
     *
     * @return parsed value of type T
     */
    template <typename T>
    T parse(const std::string& str)
    {
      using namespace std::string_literals;

      T ret;

      auto ss = std::stringstream(str);
      ss >> ret;

      if (ss.fail() || !ss.eof())
      {
        throw std::runtime_error(
            "Failed to parse \""s + str + "\" as value of type \'" + typeid(T).name() + "\'");
      }
      return ret;
    }

  } // namespace Utilities
} // namespace GridKit
