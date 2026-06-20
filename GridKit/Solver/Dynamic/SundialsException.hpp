#pragma once

#include <exception>

namespace AnalysisManager
{
  namespace Sundials
  {
    class SundialsException : public std::exception
    {
    public:
      const char* what() const noexcept override
      {
        return "Method in a SUNDIALS solver wrapper failed!\n";
      }
    };
  } // namespace Sundials
} // namespace AnalysisManager
