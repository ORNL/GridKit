#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelayData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalHistory.hpp>
#include <GridKit/Model/PhasorDynamics/SignalPorts.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit::PhasorDynamics
{
  enum class SignalDelayInternalVariables : size_t
  {
    Z, ///< Delayed output.
  };

  enum class SignalDelayExternalVariables : size_t
  {
    U, ///< Current source value.
  };

  /// Algebraic delayed signal, created by the system composer.
  template <typename scalar_type, typename index_type>
  class SignalDelay : public Component<scalar_type, index_type>
  {
    using Component<scalar_type, index_type>::gridkit_component_id_;
    using Component<scalar_type, index_type>::size_;
    using Component<scalar_type, index_type>::allocated_;
    using Component<scalar_type, index_type>::time_;
    using Component<scalar_type, index_type>::y_;
    using Component<scalar_type, index_type>::yp_;
    using Component<scalar_type, index_type>::f_;
    using Component<scalar_type, index_type>::ws_;
    using Component<scalar_type, index_type>::ws_indices_;
    using Component<scalar_type, index_type>::tag_;
    using Component<scalar_type, index_type>::abs_tol_;
    using Component<scalar_type, index_type>::variable_indices_;
    using Component<scalar_type, index_type>::residual_indices_;
    using Component<scalar_type, index_type>::nnz_;
    using Component<scalar_type, index_type>::J_rows_buffer_;
    using Component<scalar_type, index_type>::J_cols_buffer_;
    using Component<scalar_type, index_type>::J_vals_buffer_;

  public:
    using ScalarT      = scalar_type;
    using IdxT         = index_type;
    using RealT        = typename Component<ScalarT, IdxT>::RealT;
    using MonitorT     = Model::VariableMonitor<SignalDelay, SignalDelayData>;
    using ModelDataT   = SignalDelayData<RealT, IdxT>;
    using SignalPortsT = SignalPorts<ScalarT, ModelDataT>;

    SignalDelay(const ModelDataT& data, SignalHistory<RealT>& history);
    ~SignalDelay() override;
    int  setGridKitComponentID(IdxT id) override;
    int  allocate() override;
    int  verify() const override;
    int  initialize() override;
    int  tagDifferentiable() override;
    int  setAbsoluteTolerance(RealT tolerance) override;
    int  evaluateResidual() override;
    int  evaluateJacobian() override;
    void linkOutput();

    SignalPortsT& getPorts()
    {
      return ports_;
    }

    const Model::VariableMonitorBase* getMonitor() const override;

    int evaluateInternalResidual(const ScalarT* y, const ScalarT* yp, const ScalarT* buses, const ScalarT* signals, ScalarT* residual);

  private:
    RealT                     history_value_{0};
    RealT                     source_weight_{0};
    RealT                     delay_;
    SignalHistory<RealT>&     history_;
    SignalPortsT              ports_;
    std::unique_ptr<MonitorT> monitor_;
  };
} // namespace GridKit::PhasorDynamics
