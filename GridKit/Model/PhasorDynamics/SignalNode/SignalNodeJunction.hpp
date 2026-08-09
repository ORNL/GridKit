/**
 * @file SignalNodeJunction.hpp
 * @brief Internal algebraic component for a signal-node junction.
 */
#pragma once

#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Owns the algebraic output variable of a signal-node junction.
     *
     * The component enforces
     * \f$f = y - b - \sum_i g_i u_i = 0\f$. It is created internally from a
     * junction definition on a signal node and is not a user-facing device.
     */
    template <typename scalar_type, typename index_type>
    class SignalNodeJunction : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::allocated_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::y_ext_;

    public:
      using ScalarT       = scalar_type;
      using IdxT          = index_type;
      using RealT         = typename Component<ScalarT, IdxT>::RealT;
      using SignalT       = SignalNode<ScalarT, IdxT>;
      using JunctionInput = typename SignalT::JunctionInput;

      /**
       * @brief Construct and configure a junction output.
       *
       * @param[in] output Junction output signal node.
       * @param[in] bias Constant junction offset.
       * @param[in] initialization_input_index Resolved position in @p inputs
       *            to back-solve during initialization.
       * @param[in] inputs Resolved, weighted input signal nodes.
       */
      SignalNodeJunction(SignalT*                   output,
                         RealT                      bias,
                         IdxT                       initialization_input_index,
                         std::vector<JunctionInput> inputs);

      ~SignalNodeJunction() override = default;

      int setGridKitComponentID(IdxT component_id) override final;
      int allocate() override final;
      int verify() const override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int setAbsoluteTolerance(RealT rel_tol) override final;
      int evaluateInternalResidual() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

    private:
      void gatherExternalVariables();

      SignalT* output_{nullptr}; ///< Non-owning junction output
    };

  } // namespace PhasorDynamics
} // namespace GridKit
