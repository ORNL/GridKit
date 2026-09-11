/**
 * @file Port.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief A connection to a signal node.
 */

#pragma once

#include <cassert>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Representation of a connection to a @ref SignalNode.
    template <typename scalar_type, typename index_type>
    class Port
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using SignalNodeT = SignalNode<ScalarT, IdxT>;

      /// Connect a signal node to this port.
      void connect(SignalNodeT* node) noexcept
      {
        assert(node != nullptr);
        signal_node_ = node;
        assign(node);
      }

      /// Overload indicating if a signal node is connected to this port.
      operator bool() const noexcept
      {
        return connected();
      }

      /// Whether or not a signal node is connected to this port.
      bool connected() const noexcept
      {
        return signal_node_ != nullptr;
      }

      /// Whether or not the connected signal node has been linked to an output
      /// port.
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
  } // namespace PhasorDynamics
} // namespace GridKit
