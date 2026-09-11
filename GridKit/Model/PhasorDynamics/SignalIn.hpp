/**
 * @file SignalIn.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief Receiver of a signal from a signal node.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Port.hpp>

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
