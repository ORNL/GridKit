/**
 * @file ConfigurationChecks.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Collects the configuration problems found by a model.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Collects the configuration problems found by a model.
     *
     * A model's verify() fills one instance with a message for every
     * condition that does not hold and returns it. The caller decides how
     * the messages are reported; passed() is true when there are none.
     */
    class ConfigurationChecks
    {
    public:
      /// Record a problem.
      void fail(std::string message)
      {
        errors_.push_back(std::move(message));
      }

      /// Record the message when the condition does not hold.
      void check(bool condition, std::string message)
      {
        if (!condition)
        {
          fail(std::move(message));
        }
      }

      /// True when no check has failed.
      bool passed() const
      {
        return errors_.empty();
      }

      /// Every recorded message, in the order the checks ran.
      const std::vector<std::string>& errors() const
      {
        return errors_;
      }

    private:
      std::vector<std::string> errors_;
    };
  } // namespace Model
} // namespace GridKit
