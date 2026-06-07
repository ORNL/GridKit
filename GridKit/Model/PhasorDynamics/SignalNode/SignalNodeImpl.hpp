/**
 * @file SignalNode model implementation.
 */
#include <utility>

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
    void SignalNode<scalar_type, index_type>::set(ScalarT* signal, IdxT* variable_index, PrehistoryFn prehistory)
    {
      signal_         = signal;
      variable_index_ = variable_index;
      prehistory_     = std::move(prehistory);
      if (signal_)
      {
        initial_value_ = realValue();
      }
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
    scalar_type SignalNode<scalar_type, index_type>::readWithDelay(RealT delay) const
    {
      if (delay <= 0.0)
      {
        return read();
      }
      return ScalarT{history_.read(read_time_ - delay, prehistory_, initial_value_)};
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::init(ScalarT signal)
    {
      *signal_ = signal;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::requireHistoryWindow(RealT delay)
    {
      if (delay > history_window_)
      {
        history_window_ = delay;
        history_.requireWindow(delay);
      }
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::setReadTime(RealT t)
    {
      read_time_ = t;
    }

    template <typename scalar_type, typename index_type>
    void SignalNode<scalar_type, index_type>::resetHistory(RealT t0)
    {
      history_.reset();
      if (linked())
      {
        initial_value_ = realValue();
        history_.record(t0, initial_value_);
      }
    }

    template <typename scalar_type, typename index_type>
    int SignalNode<scalar_type, index_type>::stepAccepted(RealT t)
    {
      if (linked())
      {
        history_.record(t, realValue());
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    auto SignalNode<scalar_type, index_type>::realValue() const -> RealT
    {
      if constexpr (std::is_arithmetic_v<ScalarT>)
      {
        return static_cast<RealT>(*signal_);
      }
      else
      {
        return static_cast<RealT>(signal_->getValue());
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
