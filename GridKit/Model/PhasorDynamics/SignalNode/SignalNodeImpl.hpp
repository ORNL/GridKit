/**
 * @file SignalNode model implementation.
 */
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    SignalNode<scalar_type, index_type>::SignalNode()
    {
    }

    template <typename scalar_type, typename index_type>
    SignalNode<scalar_type, index_type>::SignalNode(const SignalNodeData<RealT, IdxT>& data)
      : signal_id_(data.signal_id)
    {
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::set(ScalarT* signal, IdxT* variable_index)
    {
      signal_         = signal;
      variable_index_ = variable_index;
    }

    template <typename scalar_type, typename index_type>
    bool SignalNode<scalar_type, index_type>::linked() const
    {
      return (signal_) && (variable_index_);
    }

    template <typename scalar_type, typename index_type>
    scalar_type SignalNode<scalar_type, index_type>::read() const
    {
      return *signal_;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::init(ScalarT signal)
    {
      *signal_ = signal;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
