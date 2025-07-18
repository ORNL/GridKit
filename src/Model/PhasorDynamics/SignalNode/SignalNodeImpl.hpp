/**
 * @file SignalNode model implementation.
 */
#include "SignalNode.hpp"
#include <Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    SignalNode<ScalarT, IdxT>::SignalNode()
      : size_(0)
    {
    }

    template <class ScalarT, typename IdxT>
    SignalNode<ScalarT, IdxT>::SignalNode(const SignalNodeData<real_type, IdxT>& data)
      : signal_id_(data.signal_id),
        size_(0)
    {
    }

    template <class ScalarT, typename IdxT>
    int SignalNode<ScalarT, IdxT>::allocate()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int SignalNode<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int SignalNode<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int SignalNode<ScalarT, IdxT>::evaluateResidual()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int SignalNode<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    void SignalNode<ScalarT, IdxT>::set(ScalarT* signal)
    {
      signal_ = signal;
    }

    template <class ScalarT, typename IdxT>
    ScalarT SignalNode<ScalarT, IdxT>::read() const
    {
      return *signal_;
    }

    template <class ScalarT, typename IdxT>
    void SignalNode<ScalarT, IdxT>::init(ScalarT signal)
    {
      *signal_ = signal;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
