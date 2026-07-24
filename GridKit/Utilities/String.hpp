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

    std::string strip(std::string str)
    {
      auto notspace = [](char c)
      { return !std::isspace(c); };
      str.erase(begin(str), std::find_if(begin(str), end(str), notspace));
      str.erase(std::find_if(rbegin(str), rend(str), notspace).base(), end(str));
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

      auto ss = std::stringstream(strip(str));
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
