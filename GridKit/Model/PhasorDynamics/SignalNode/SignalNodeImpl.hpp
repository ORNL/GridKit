/**
 * @file SignalNode model implementation.
 */
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    SignalNode<ScalarP, IdxP>::SignalNode()
    {
    }

    template <typename ScalarP, typename IdxP>
    SignalNode<ScalarP, IdxP>::SignalNode(const SignalNodeData<RealT, IdxT>& data)
      : signal_id_(data.signal_id)
    {
    }

    template <typename ScalarP, typename IdxP>
    void SignalNode<ScalarP, IdxP>::set(ScalarT* signal, IdxT* variable_index)
    {
      signal_         = signal;
      variable_index_ = variable_index;
    }

    template <typename ScalarP, typename IdxP>
    bool SignalNode<ScalarP, IdxP>::linked() const
    {
      return (signal_) && (variable_index_);
    }

    template <typename ScalarP, typename IdxP>
    ScalarP SignalNode<ScalarP, IdxP>::read() const
    {
      return *signal_;
    }

    template <typename ScalarP, typename IdxP>
    void SignalNode<ScalarP, IdxP>::init(ScalarT signal)
    {
      *signal_ = signal;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
