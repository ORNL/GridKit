/**
 * @file SignalNode model implementation.
 */

#include <cassert>

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
    index_type SignalNode<scalar_type, index_type>::signalId() const noexcept
    {
      return signal_id_;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::setAssigned() noexcept
    {
      assert(!assigned_);
      assigned_ = true;
    }

    template <typename scalar_type, typename index_type>
    bool SignalNode<scalar_type, index_type>::assigned() const noexcept
    {
      return assigned_;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::link(ScalarT* signal, IdxT* variable_index) noexcept
    {
      signal_         = signal;
      variable_index_ = variable_index;
    }

    template <typename scalar_type, typename index_type>
    bool SignalNode<scalar_type, index_type>::linked() const noexcept
    {
      return (signal_) && (variable_index_);
    }

    template <typename scalar_type, typename index_type>
    scalar_type SignalNode<scalar_type, index_type>::read() const noexcept
    {
      assert(signal_);
      return *signal_;
    }

    template <typename scalar_type, typename index_type>
    index_type SignalNode<scalar_type, index_type>::getVariableIndex() const noexcept
    {
      assert(variable_index_);
      return *variable_index_;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::init(ScalarT signal) noexcept
    {
      assert(signal_);
      *signal_ = signal;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
