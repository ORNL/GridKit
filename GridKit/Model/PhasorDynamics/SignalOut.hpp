#pragma once

#include <GridKit/Model/PhasorDynamics/Port.hpp>

namespace GridKit::PhasorDynamics
{
  /// Port for sending a signal, making it available via a @ref SignalNode.
  template <typename scalar_type, typename index_type>
  class SignalOut : public Port<scalar_type, index_type>
  {
  public:
    using ScalarT     = scalar_type;
    using IdxT        = index_type;
    using SignalNodeT = SignalNode<ScalarT, IdxT>;

    void link(ScalarT* signal_var, IdxT* global_index)
    {
      assert(this->connected());
      this->signal_node_->link(signal_var, global_index);
    }

  protected:
    void assign(SignalNodeT* node) const override
    {
      node->setAssigned();
    }
  };
} // namespace GridKit::PhasorDynamics
