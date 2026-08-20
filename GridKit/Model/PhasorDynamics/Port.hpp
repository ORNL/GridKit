#pragma once

#include <cassert>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>

namespace GridKit::PhasorDynamics
{
  /// Representation of a connection to a @ref SignalNode.
  template <typename scalar_type, typename index_type>
  class Port
  {
  public:
    using ScalarT     = scalar_type;
    using IdxT        = index_type;
    using SignalNodeT = SignalNode<ScalarT, IdxT>;

    void connect(SignalNodeT* node) noexcept
    {
      assert(node != nullptr);
      signal_node_ = node;
      assign(node);
    }

    operator bool() const noexcept
    {
      return connected();
    }

    bool connected() const noexcept
    {
      return signal_node_ != nullptr;
    }

    bool linked() const noexcept
    {
      return connected() ? signal_node_->linked() : false;
    }

  protected:
    virtual void assign([[maybe_unused]] SignalNodeT*) const
    {
    }

    SignalNodeT* signal_node_{nullptr};
  };
} // namespace GridKit::PhasorDynamics
