#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename real_type, typename index_type>
    struct SignalNodeData;
  } // namespace PhasorDynamics
} // namespace GridKit

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
     * (Component):[OutputPort] -> {SignalNode} -> [InputPort]:(Component)
     *
     * A SignalNode can be "connected" to a Port. When that port is an
     * OutputPort, the SignalNode is considered `assigned()` since it can be
     * connected to only one OutputPort. The SignalNode is considered `linked()`
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

      SignalNode();
      SignalNode(const SignalNodeData<RealT, IdxT>& data);

      virtual ~SignalNode() = default;

      IdxT    signalId() const noexcept;
      void    setAssigned() noexcept;
      bool    assigned() const noexcept;
      void    link(ScalarT* signal_in, IdxT* global_index) noexcept;
      bool    linked() const noexcept;
      ScalarT read() const noexcept;
      IdxT    getVariableIndex() const noexcept;
      void    init(ScalarT signal_in) noexcept;

    private:
      ScalarT* signal_{nullptr};
      IdxT     signal_id_{0};
      IdxT*    variable_index_{nullptr};
      bool     assigned_{false};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
