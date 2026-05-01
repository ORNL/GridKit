#pragma once

#include <algorithm>
#include <cctype>
#include <string>

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Convert a string to all uppercase
     */
    std::string toUpper(std::string str)
    {
      std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c)
                     { return std::toupper(c); });
      return str;
    }

  } // namespace Utilities
} // namespace GridKit
