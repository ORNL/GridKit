#pragma once

#include <cassert>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeData.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief SignalNode model implementation base class.
     *
     * A signal node conceptually carries a signal (scalar value) from an
     * output port of one component to an input port of another component (or
     * multiple components).
     *
     * (Component):[SignalOut] -> {SignalNode} -> [SignalIn]:(Component)
     *
     * A SignalNode can be "connected" to a Port. When that port is an
     * SignalOut, the SignalNode is considered `assigned()` since it can be
     * connected to only one SignalOut. The SignalNode is considered `linked()`
     * when the actual signal (scalar variable) has been made available.
     *
     */
    template <typename scalar_type, typename index_type>
    class SignalNode
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      SignalNode()
      {
      }

      SignalNode(const SignalNodeData<RealT, IdxT>& data)
        : signal_id_(data.signal_id)
      {
      }

      virtual ~SignalNode() = default;

      [[gnu::always_inline]]
      IdxT signalId() const noexcept
      {
        return signal_id_;
      }

      [[gnu::always_inline]]
      void setAssigned() noexcept
      {
        assert(!assigned_);
        assigned_ = true;
      }

      [[gnu::always_inline]]
      bool assigned() const noexcept
      {
        return assigned_;
      }

      [[gnu::always_inline]]
      void link(ScalarT* signal_in, IdxT* global_index) noexcept
      {
        signal_         = signal_in;
        variable_index_ = global_index;
      }

      [[gnu::always_inline]]
      bool linked() const noexcept
      {
        return (signal_) && (variable_index_);
      }

      [[gnu::always_inline]]
      ScalarT read() const noexcept
      {
        assert(signal_);
        return *signal_;
      }

      [[gnu::always_inline]]
      IdxT getVariableIndex() const noexcept
      {
        assert(variable_index_);
        return *variable_index_;
      }

      [[gnu::always_inline]]
      void init(ScalarT signal_in) noexcept
      {
        assert(signal_);
        *signal_ = signal_in;
      }

    private:
      ScalarT* signal_{nullptr};
      IdxT     signal_id_{0};
      IdxT*    variable_index_{nullptr};
      bool     assigned_{false};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
