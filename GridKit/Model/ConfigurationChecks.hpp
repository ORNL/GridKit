/**
 * @file ConfigurationChecks.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Collects the configuration problems found by a model.
 */

#pragma once

#include <format>
#include <string>
#include <string_view>
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
    class [[nodiscard]] ConfigurationChecks
    {
    public:
      ConfigurationChecks() = default;

      /// @param[in] model Model name that prefixes every recorded message.
      explicit ConfigurationChecks(std::string_view model)
        : prefix_(std::format("{}: ", model))
      {
      }

      /// Record a problem.
      void fail(std::string_view message)
      {
        errors_.push_back(prefix_ + std::string(message));
      }

      /// Record the message when the condition does not hold.
      void check(bool condition, std::string_view message)
      {
        if (!condition)
        {
          fail(message);
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
      std::string              prefix_;
      std::vector<std::string> errors_;
    };
  } // namespace Model
} // namespace GridKit
