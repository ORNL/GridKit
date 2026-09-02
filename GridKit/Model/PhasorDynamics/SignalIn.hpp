/**
 * @file SignalIn.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief Receiver of a signal from a signal node.
 */

#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Port.hpp>
#include <GridKit/Utilities/ConfigurationChecks.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Port for receiving a signal from a @ref SignalNode.
    template <typename scalar_type, typename index_type>
    class SignalIn : public Port<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

      /// Read a value from the connected signal node.
      ScalarT readSignal() const
      {
        assert(this->connected());
        return this->signal_node_->read();
      }

      /// Retrieve the variable index of the connected signal node.
      IdxT signalVariableIndex() const
      {
        assert(this->connected());
        return this->signal_node_->getVariableIndex();
      }

      /// Verify an optional input is linked whenever it is connected.
      void checkOptional(Utilities::ConfigurationChecks& checks, const char* name) const
      {
        if (this->connected() && !this->linked())
        {
          checks.fail() << name << " signal attached with no linked source\n";
        }
      }

      /// Verify a required input is connected and linked.
      void checkRequired(Utilities::ConfigurationChecks& checks, const char* name) const
      {
        if (!this->connected())
        {
          checks.fail() << name << " signal is required\n";
          return;
        }
        checkOptional(checks, name);
      }

      /// Read the connected signal, or the fallback when disconnected.
      ScalarT readOrDefault(ScalarT fallback) const
      {
        if (this->connected())
        {
          return readSignal();
        }
        return fallback;
      }

      /// Refresh an explicitly selected workspace value and global index.
      void refreshWorkspace(ScalarT fallback, ScalarT& value, IdxT& index) const
      {
        value = fallback;
        index = INVALID_INDEX<IdxT>;
        if (this->connected())
        {
          value = readSignal();
          index = signalVariableIndex();
        }
      }

      /// Write a value to the connected signal node.
      ///
      /// @warning Use only during initialization as this violates assumptions.
      void writeValue(ScalarT value)
      {
        this->signal_node_->init(value);
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit
