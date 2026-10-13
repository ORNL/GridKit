#pragma once

#include <GridKit/Model/PhasorDynamics/SignalDelay/SignalDelay.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit::PhasorDynamics
{
  template <typename scalar_type, typename index_type>
  SignalDelay<scalar_type, index_type>::~SignalDelay() = default;

  template <typename scalar_type, typename index_type>
  const Model::VariableMonitorBase* SignalDelay<scalar_type, index_type>::getMonitor() const
  {
    return monitor_.get();
  }

  template <typename scalar_type, typename index_type>
  SignalDelay<scalar_type, index_type>::SignalDelay(const ModelDataT& data, SignalHistory<RealT>& history)
    : delay_(std::get<RealT>(data.parameters.at(SignalDelayParameters::tau))), history_(history)
  {
    if (!std::isfinite(delay_) || delay_ <= 0)
    {
      throw std::invalid_argument("SignalDelay requires a finite positive delay");
    }
    size_    = 1;
    monitor_ = std::make_unique<MonitorT>(data);
    monitor_->set(SignalDelayMonitors::z, [this]
                  { return y_.getData()[static_cast<size_t>(SignalDelayInternalVariables::Z)]; });
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::setGridKitComponentID(IdxT id)
  {
    gridkit_component_id_ = id;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  void SignalDelay<scalar_type, index_type>::linkOutput()
  {
    ports_.out.template port<SignalDelayOutputs::z>().link(y_.getData(), &this->getVariableIndex(static_cast<IdxT>(SignalDelayInternalVariables::Z)));
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::allocate()
  {
    if (!allocated_)
    {
      this->allocateVectors(1);
    }
    ws_.resize(1);
    ws_indices_.assign(1, INVALID_INDEX<IdxT>);
    variable_indices_.assign(1, 0);
    residual_indices_.assign(1, 0);
    tag_.assign(1, false);
    linkOutput();
    allocated_ = true;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::verify() const
  {
    if (!ports_.in.template port<SignalDelayInputs::u>().linked() || !ports_.out.template port<SignalDelayOutputs::z>().linked())
    {
      throw std::invalid_argument("SignalDelay requires linked input and output signals");
    }
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::initialize()
  {
    const auto Z             = static_cast<size_t>(SignalDelayInternalVariables::Z);
    const auto interpolation = history_.sample(time_, delay_);
    y_.getData()[Z]          = interpolation.constant
                      + interpolation.weight * ports_.in.template port<SignalDelayInputs::u>().readSignal();
    yp_.getData()[Z] = history_.initialDerivative(time_ - delay_);
    y_.setDataUpdated();
    yp_.setDataUpdated();
    linkOutput();
    if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
    {
      this->initializeDependencyTrackingVariableNumbers();
    }
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::tagDifferentiable()
  {
    tag_[static_cast<size_t>(SignalDelayInternalVariables::Z)] = false;
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
  {
    abs_tol_.getData()[static_cast<size_t>(SignalDelayInternalVariables::Z)] = tolerance;
    abs_tol_.setDataUpdated();
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::evaluateResidual()
  {
    const auto U     = static_cast<size_t>(SignalDelayExternalVariables::U);
    const auto input = ports_.in.template port<SignalDelayInputs::u>();
    auto*      ws    = ws_.getData();
    ws[U]            = input.readSignal();
    ws_indices_[U]   = input.signalVariableIndex();
    ws_.setDataUpdated();

    const auto interpolation = history_.sample(time_, delay_);
    history_value_           = interpolation.constant;
    source_weight_           = interpolation.weight;
    evaluateInternalResidual(y_.getData(), yp_.getData(), nullptr, ws, f_.getData());
    f_.setDataUpdated();
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::evaluateInternalResidual(
      const ScalarT* y, const ScalarT*, const ScalarT*, const ScalarT* signals, ScalarT* residual)
  {
    const auto Z = static_cast<size_t>(SignalDelayInternalVariables::Z);
    const auto U = static_cast<size_t>(SignalDelayExternalVariables::U);
    residual[Z]  = -y[Z] + history_value_ + source_weight_ * signals[U];
    return 0;
  }

  template <typename scalar_type, typename index_type>
  int SignalDelay<scalar_type, index_type>::evaluateJacobian()
  {
    if (J_rows_buffer_ == nullptr)
    {
      nnz_           = 2;
      J_rows_buffer_ = new IdxT[2];
      J_cols_buffer_ = new IdxT[2];
      J_vals_buffer_ = new RealT[2];
    }
    const auto interpolation = history_.sample(time_, delay_);
    const auto source        = ports_.in.template port<SignalDelayInputs::u>().signalVariableIndex();
    const auto Z             = static_cast<IdxT>(SignalDelayInternalVariables::Z);
    J_rows_buffer_[0] = J_rows_buffer_[1] = this->getResidualIndex(Z);
    J_cols_buffer_[0]                     = this->getVariableIndex(Z);
    J_vals_buffer_[0]                     = -1;
    // Preserve the structural entry when the query moves into accepted history.
    J_cols_buffer_[1]                     = source == INVALID_INDEX<IdxT> ? this->getVariableIndex(Z) : source;
    J_vals_buffer_[1]                     = source == INVALID_INDEX<IdxT> ? 0 : interpolation.weight;
    return this->constructCoo();
  }
} // namespace GridKit::PhasorDynamics
