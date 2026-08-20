#pragma once

#include <GridKit/Model/PhasorDynamics/Port.hpp>

namespace GridKit::PhasorDynamics
{
  /// Port for receiving a signal from a @ref SignalNode.
  template <typename scalar_type, typename index_type>
  class SignalIn : public Port<scalar_type, index_type>
  {
  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;

    ScalarT readSignal() const
    {
      assert(this->connected());
      return this->signal_node_->read();
    }

    IdxT signalVariableIndex() const
    {
      assert(this->connected());
      return this->signal_node_->getVariableIndex();
    }

    /// Use only during initialization.
    void writeValue(ScalarT value)
    {
      this->signal_node_->init(value);
    }
  };
} // namespace GridKit::PhasorDynamics
