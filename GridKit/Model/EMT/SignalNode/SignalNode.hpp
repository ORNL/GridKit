#pragma once

#include <array>

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct SignalNodeData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /*!
     * @brief SignalNode model implementation base class.
     *
     * A fully bound node exposes the value, derivative, and residual row of
     * one variable owned by another component, along with the global variable
     * and residual indices. Consumers read the value and derivative through
     * the node and accumulate external residual contributions into the
     * owner's residual row.
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

      void set(ScalarT* signal, IdxT* variable_index);
      void set(ScalarT* signal,
               ScalarT* derivative,
               ScalarT* residual,
               IdxT*    variable_index,
               IdxT*    residual_index);

      bool linked() const;
      bool derivativeLinked() const;
      bool residualLinked() const;

      ScalarT read() const;
      ScalarT readDerivative() const;
      void    accumulateResidual(ScalarT value);
      void    init(ScalarT signal_in);
      void    initDerivative(ScalarT derivative_in);

      void markDerivativeCoupling();
      bool hasDerivativeCoupling() const;

      const IdxT signalId() const
      {
        return signal_id_;
      }

      IdxT getVariableIndex() const
      {
        return *variable_index_;
      }

      IdxT getResidualIndex() const
      {
        return *residual_index_;
      }

    private:
      ScalarT* signal_{nullptr};
      ScalarT* derivative_{nullptr};
      ScalarT* residual_{nullptr};
      IdxT     signal_id_{0};

    protected:
      IdxT* variable_index_{nullptr};
      IdxT* residual_index_{nullptr};
      bool  derivative_coupling_{false};
    };

    /*!
     * @brief Three-phase electrical port bundling one signal node per phase.
     *
     * Owned by value by the component that owns the phase variables and
     * residual rows.
     */
    template <typename scalar_type, typename index_type>
    struct Port3
    {
      using NodeT = SignalNode<scalar_type, index_type>;

      std::array<NodeT, 3> nodes{};

      NodeT* a()
      {
        return &nodes[0];
      }

      NodeT* b()
      {
        return &nodes[1];
      }

      NodeT* c()
      {
        return &nodes[2];
      }
    };

  } // namespace EMT
} // namespace GridKit
