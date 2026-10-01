#pragma once

#include <stdexcept>
#include <string>

namespace GridKit
{
  namespace Utilities
  {
    class NotImplementedError : public std::runtime_error
    {
    public:
      NotImplementedError(const std::string& funcname)
        : std::runtime_error("ERROR: Function not implemented: " + funcname)
      {
      }
    };
  } // namespace Utilities
} // namespace GridKit
