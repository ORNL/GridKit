/**
 * @file SignalOut.hpp
 * @author superwhiskers <whiskerdev@protonmail.com>
 * @brief Provider of a signal to a signal node.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Port.hpp>
#include <GridKit/Model/PhasorDynamics/SignalIn.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Port for sending a signal, making it available via a @ref SignalNode.
    template <typename scalar_type, typename index_type>
    class SignalOut : public Port<scalar_type, index_type>
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using SignalNodeT = SignalNode<ScalarT, IdxT>;

      /// Link a signal variable, providing its index to the signal node. A bus
      /// voltage also provides the residual index that receives injections.
      void link(ScalarT* signal_var, IdxT* global_index, IdxT* residual_index = nullptr)
      {
        assert(this->connected());
        this->signal_node_->link(signal_var, global_index, residual_index);
      }

      /// Link the signal read by `source`, with a residual index of the caller.
      void link(const SignalIn<ScalarT, IdxT>& source, IdxT* residual_index)
      {
        assert(this->connected() && source.linked());
        this->signal_node_->link(*source.node(), residual_index);
      }

    protected:
      void assign(SignalNodeT* node) const override
      {
        node->setAssigned();
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit
