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

      /// Read the value of the signal node if connected.
      const ScalarT& readSignal() const
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

      /// Retrieve the residual index of the connected signal node.
      IdxT signalResidualIndex() const
      {
        assert(this->connected());
        return this->signal_node_->getResidualIndex();
      }

      /// Copy the signal value and its variable index.
      void gather(ScalarT& value, IdxT& variable_index) const
      {
        value          = readSignal();
        variable_index = signalVariableIndex();
      }

      /// Copy the signal value, its variable index, and its residual index.
      void gather(ScalarT& value, IdxT& variable_index, IdxT& residual_index) const
      {
        gather(value, variable_index);
        residual_index = signalResidualIndex();
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
