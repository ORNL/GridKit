#pragma once

namespace GridKit
{
  namespace Utilities
  {
    namespace Colors
    {
      // must be const pointer and const dest for
      // const string declarations to pass -Wwrite-strings
      static const char* const RED    = "\033[1;31m";
      static const char* const GREEN  = "\033[1;32m";
      static const char* const YELLOW = "\033[33;1m";
      static const char* const BLUE   = "\033[34;1m";
      static const char* const ORANGE = "\u001b[38;5;208m";
      static const char* const CLEAR  = "\033[0m";
    } // namespace Colors
  } // namespace Utilities
} // namespace GridKit
