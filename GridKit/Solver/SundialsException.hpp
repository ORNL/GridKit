#pragma once

#include <exception>
#include <string>

namespace AnalysisManager
{
  namespace Sundials
  {
    /// Exception wrapping an error that occurs within SUNDIALS.
    class SundialsException : public std::exception
    {
      std::string message_;

    public:
      explicit SundialsException(
        const char* message = "A SUNDIALS class wrapper failed!\n"
      )
        : message_(message)
      {
      }

      const char* what() const noexcept override
      {
        return message_.c_str();
      }
    };
  } // namespace Sundials
} // namespace AnalysisManager
