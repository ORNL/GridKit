/**
 * @file ConfigurationChecks.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Error accumulator for model configuration validation.
 */

#pragma once

#include <ostream>

#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Utilities
  {
    /**
     * @brief Error accumulator for model configuration validation.
     *
     * One instance collects the errors found while loading parameters or
     * verifying a model, logging each with the model-name prefix. The
     * caller reports errorCount() as its error total.
     */
    class ConfigurationChecks
    {
    public:
      explicit ConfigurationChecks(const char* model)
        : model_(model)
      {
      }

      /// Log one error against this model and count it.
      std::ostream& fail()
      {
        ++error_count_;
        return Logger::error() << model_ << ": ";
      }

      /// Log and count one error when the condition does not hold.
      void check(bool condition, const char* message)
      {
        if (!condition)
        {
          fail() << message << '\n';
        }
      }

      int errorCount() const
      {
        return error_count_;
      }

    private:
      const char* model_;
      int         error_count_{0};
    };
  } // namespace Utilities
} // namespace GridKit
