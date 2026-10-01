#pragma once

#include <iomanip>
#include <iostream>

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Print `width` number of spaces
     */
    void leftPad(std::ostream& os, unsigned width)
    {
      os << std::setfill(' ') << std::setw(static_cast<int>(width)) << "";
    }

  } // namespace Utilities
} // namespace GridKit
