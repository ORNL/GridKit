/**
 * @file Port.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief A connection to a signal node or a bus.
 */

#pragma once

#include <cassert>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Representation of a connection to a @ref SignalNode or a bus.
    template <typename scalar_type,
              typename index_type,
              typename node_type = SignalNode<scalar_type, index_type>>
    class Port
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using NodeT   = node_type;

      /// Connect a node to this port.
      void connect(NodeT* node) noexcept
      {
        assert(node != nullptr);
        node_ = node;
        assign(node);
      }

      /// Overload indicating if a node is connected to this port.
      operator bool() const noexcept
      {
        return connected();
      }

      /// Whether or not a node is connected to this port.
      bool connected() const noexcept
      {
        return node_ != nullptr;
      }

      /// Whether or not the connected signal node has been linked to an output
      /// port.
      bool linked() const noexcept
      {
        return connected() ? node_->linked() : false;
      }

      /// The node connected to this port.
      NodeT* node() const noexcept
      {
        assert(connected());
        return node_;
      }

    protected:
      virtual void assign([[maybe_unused]] NodeT*) const
      {
      }

      NodeT* node_{nullptr};
    };
  } // namespace PhasorDynamics
} // namespace GridKit
